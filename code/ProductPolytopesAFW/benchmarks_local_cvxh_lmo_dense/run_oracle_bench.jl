#!/usr/bin/env julia
using CSV
using DataFrames
using FrankWolfe
using LinearAlgebra
using ProductPolytopesAFW
using Random
using YAML

const VALID_BACKENDS = Set([
    "vec_best",
    "mat_opt_cacheauto",
    "mat_opt_cacheoff",
    "mat_scan_cacheauto",
    "mat_scan_cacheoff",
])

function usage()
    return "Usage: julia --project=. benchmarks_local_cvxh_lmo_dense/run_oracle_bench.jl <config.yml> <n> <backend_tag>"
end

function load_config(path::AbstractString)
    isfile(path) || error("Config file not found: $path")
    cfg = YAML.load_file(path)
    get(cfg, "k", nothing) == 2 || error("This harness currently supports only k = 2.")
    return cfg
end

function parse_args(args::Vector{String})
    length(args) == 3 || error(usage())
    cfg_path = abspath(args[1])
    n = parse(Int, args[2])
    backend_tag = args[3]
    n > 0 || error("n must be positive.")
    backend_tag in VALID_BACKENDS || error("Unknown backend tag '$backend_tag'.")
    return cfg_path, n, backend_tag
end

function results_dir(cfg_path::AbstractString, cfg)
    return abspath(joinpath(dirname(cfg_path), String(cfg["results_dir"])))
end

function case_id(n::Int, backend_tag::AbstractString)
    return "oracle_n$(n)_$(backend_tag)"
end

function backend_metadata(tag::AbstractString)
    if tag == "vec_best"
        return (backend="vector", optimized="na", cache="na", cache_cap=nothing, use_optimized=false)
    elseif tag == "mat_opt_cacheauto"
        return (backend="matrix", optimized="true", cache="on", cache_cap=nothing, use_optimized=true)
    elseif tag == "mat_opt_cacheoff"
        return (backend="matrix", optimized="true", cache="off", cache_cap=0, use_optimized=true)
    elseif tag == "mat_scan_cacheauto"
        return (backend="matrix", optimized="false", cache="on", cache_cap=nothing, use_optimized=false)
    elseif tag == "mat_scan_cacheoff"
        return (backend="matrix", optimized="false", cache="off", cache_cap=0, use_optimized=false)
    end
    error("Unsupported backend tag '$tag'.")
end

function sample_polytope_matrix(
    rng::AbstractRNG,
    n_points::Int,
    n::Int,
    idx::Int,
    vertex_sampling::AbstractString,
    rho::Float64,
)
    center = (8.0 * (idx - 1)) .+ 0.5 .* randn(rng, n)
    half_lengths = 0.5 .+ rand(rng, n)
    vertices = Matrix{Float64}(undef, n_points, n)

    if vertex_sampling == "box_uniform"
        @inbounds for i in 1:n_points
            vertices[i, :] .= center .+ (2 .* rand(rng, n) .- 1.0) .* half_lengths
        end
    elseif vertex_sampling == "sphere"
        radius = rho * minimum(half_lengths)
        @inbounds for i in 1:n_points
            direction = randn(rng, n)
            norm_direction = norm(direction)
            while norm_direction == 0.0
                direction = randn(rng, n)
                norm_direction = norm(direction)
            end
            vertices[i, :] .= center .+ radius .* (direction ./ norm_direction)
        end
    elseif vertex_sampling == "ellipsoid"
        axes = rho .* half_lengths
        @inbounds for i in 1:n_points
            direction = randn(rng, n)
            norm_direction = norm(direction)
            while norm_direction == 0.0
                direction = randn(rng, n)
                norm_direction = norm(direction)
            end
            vertices[i, :] .= center .+ axes .* (direction ./ norm_direction)
        end
    else
        error("Unsupported vertex_sampling '$vertex_sampling'.")
    end

    return vertices
end

function generate_point_clouds(cfg, n::Int)
    rng = MersenneTwister(Int(cfg["seed"]) + n)
    n_points = Int(cfg["n_points_multiplier"]) * n
    sampling = String(cfg["vertex_sampling"])
    rho = Float64(cfg["sphere_radius_factor"])
    return [
        sample_polytope_matrix(rng, n_points, n, idx, sampling, rho)
        for idx in 1:Int(cfg["k"])
    ]
end

function generate_direction_batches(cfg, n::Int)
    rng = MersenneTwister(Int(cfg["seed"]) + 10_000 + n)
    batch_size = Int(cfg["direction_batch_size"])
    return [[randn(rng, n) for _ in 1:batch_size] for _ in 1:Int(cfg["k"])]
end

function build_lmos(vertices::Vector{Matrix{Float64}}, backend_tag::AbstractString)
    meta = backend_metadata(backend_tag)
    if meta.backend == "vector"
        return [
            ProductPolytopesAFW.CopyExtremePointLMO(FrankWolfe.ConvexHullLMO(collect(eachrow(V))))
            for V in vertices
        ], meta
    end

    return [
        ProductPolytopesAFW.MatrixConvexHullLMO(
            V;
            cache_cap=meta.cache_cap,
            use_optimized_search=meta.use_optimized,
        )
        for V in vertices
    ], meta
end

function run_oracle_workload(vertices, directions, cfg, backend_tag::AbstractString)
    lmos, _ = build_lmos(vertices, backend_tag)
    repetitions = Int(cfg["oracle_repetitions"])
    checksum = 0.0
    @inbounds for _ in 1:repetitions
        for block_idx in eachindex(lmos)
            lmo = lmos[block_idx]
            for direction in directions[block_idx]
                atom = FrankWolfe.compute_extreme_point(lmo, direction)
                checksum += atom[1]
            end
        end
    end
    return checksum
end

function measure(workload::Function)
    GC.gc()
    alloc_bytes = @allocated workload()
    GC.gc()
    alloc_count = @allocations workload()
    GC.gc()
    value = 0.0
    elapsed_s = @elapsed value = workload()
    return value, elapsed_s, alloc_bytes, alloc_count
end

function main(args::Vector{String})
    cfg_path, n, backend_tag = parse_args(args)
    cfg = load_config(cfg_path)
    out_dir = results_dir(cfg_path, cfg)
    mkpath(joinpath(out_dir, "csv"))

    vertices = generate_point_clouds(cfg, n)
    directions = generate_direction_batches(cfg, n)
    workload = () -> run_oracle_workload(vertices, directions, cfg, backend_tag)
    checksum, elapsed_s, alloc_bytes, alloc_count = measure(workload)
    meta = backend_metadata(backend_tag)
    n_points = size(vertices[1], 1)

    row = (
        case_id=case_id(n, backend_tag),
        benchmark_kind="oracle",
        n=n,
        n_points=n_points,
        backend=meta.backend,
        optimized=meta.optimized,
        cache=meta.cache,
        line_search="na",
        iterations=Int(cfg["oracle_repetitions"]) * Int(cfg["direction_batch_size"]) * Int(cfg["k"]),
        elapsed_s=elapsed_s,
        alloc_bytes=alloc_bytes,
        alloc_count=alloc_count,
        final_primal=NaN,
        checksum=checksum,
        status="success",
        error_message="",
    )

    output_file = joinpath(out_dir, "csv", "$(row.case_id).csv")
    CSV.write(output_file, DataFrame([row]))
    println("Wrote $(output_file)")
end

try
    main(ARGS)
catch err
    showerror(stderr, err)
    println(stderr)
    println(stderr, usage())
    exit(1)
end
