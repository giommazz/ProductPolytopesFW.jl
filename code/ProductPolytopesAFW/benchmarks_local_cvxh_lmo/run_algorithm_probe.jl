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
const VALID_ALGORITHMS = Set(["fw", "afw"])

function usage()
    return "Usage: julia --project=. benchmarks_local_cvxh_lmo/run_algorithm_probe.jl <config.yml> <fw|afw> <n> <backend_tag>"
end

function load_config(path::AbstractString)
    isfile(path) || error("Config file not found: $path")
    cfg = YAML.load_file(path)
    get(cfg, "k", nothing) == 2 || error("This harness currently supports only k = 2.")
    return cfg
end

function parse_args(args::Vector{String})
    length(args) == 4 || error(usage())
    cfg_path = abspath(args[1])
    algorithm = args[2]
    n = parse(Int, args[3])
    backend_tag = args[4]
    algorithm in VALID_ALGORITHMS || error("Unknown algorithm '$algorithm'.")
    n > 0 || error("n must be positive.")
    backend_tag in VALID_BACKENDS || error("Unknown backend tag '$backend_tag'.")
    return cfg_path, algorithm, n, backend_tag
end

function results_dir(cfg_path::AbstractString, cfg)
    return abspath(joinpath(dirname(cfg_path), String(cfg["results_dir"])))
end

function case_id(algorithm::AbstractString, n::Int, backend_tag::AbstractString)
    return "$(algorithm)_n$(n)_$(backend_tag)"
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

function starting_point(lmos, n::Int)
    zero_dir = zeros(Float64, n)
    blocks = [Vector{Float64}(FrankWolfe.compute_extreme_point(lmo, zero_dir)) for lmo in lmos]
    return FrankWolfe.BlockVector(blocks)
end

function run_probe(cfg, algorithm::AbstractString, vertices, backend_tag::AbstractString)
    lmos, _ = build_lmos(vertices, backend_tag)
    prod_lmo = FrankWolfe.ProductLMO(Tuple(lmos))
    x0 = starting_point(lmos, size(vertices[1], 2))
    line_search = ProductPolytopesAFW.SafeGoldenratio(1e-9)
    iterations = algorithm == "fw" ? Int(cfg["fw_iterations"]) : Int(cfg["afw_iterations"])

    if algorithm == "fw"
        result = FrankWolfe.frank_wolfe(
            ProductPolytopesAFW.convex_feasibility_objective_v2b,
            ProductPolytopesAFW.convex_feasibility_gradient_v2!,
            prod_lmo,
            x0;
            epsilon=0.0,
            max_iteration=iterations,
            line_search=line_search,
            print_iter=typemax(Int),
            memory_mode=FrankWolfe.InplaceEmphasis(),
            verbose=false,
            trajectory=false,
        )
    else
        result = FrankWolfe.away_frank_wolfe(
            ProductPolytopesAFW.convex_feasibility_objective_v2b,
            ProductPolytopesAFW.convex_feasibility_gradient_v2!,
            prod_lmo,
            x0;
            epsilon=0.0,
            max_iteration=iterations,
            line_search=line_search,
            print_iter=typemax(Int),
            memory_mode=FrankWolfe.InplaceEmphasis(),
            verbose=false,
            trajectory=false,
        )
    end

    return result, iterations
end

function measure(workload::Function)
    GC.gc()
    alloc_bytes = @allocated workload()
    GC.gc()
    alloc_count = @allocations workload()
    GC.gc()
    value = nothing
    elapsed_s = @elapsed value = workload()
    return value, elapsed_s, alloc_bytes, alloc_count
end

function success_row(cfg, algorithm::AbstractString, n::Int, backend_tag::AbstractString, elapsed_s, alloc_bytes, alloc_count, result, iterations)
    meta = backend_metadata(backend_tag)
    n_points = Int(cfg["n_points_multiplier"]) * n
    return (
        case_id=case_id(algorithm, n, backend_tag),
        benchmark_kind=algorithm,
        n=n,
        n_points=n_points,
        backend=meta.backend,
        optimized=meta.optimized,
        cache=meta.cache,
        line_search="safegoldenratio",
        iterations=iterations,
        elapsed_s=elapsed_s,
        alloc_bytes=alloc_bytes,
        alloc_count=alloc_count,
        final_primal=Float64(result.primal),
        checksum=NaN,
        status="success",
        error_message="",
    )
end

function failure_row(cfg, algorithm::AbstractString, n::Int, backend_tag::AbstractString, message::AbstractString)
    meta = backend_metadata(backend_tag)
    n_points = Int(cfg["n_points_multiplier"]) * n
    return (
        case_id=case_id(algorithm, n, backend_tag),
        benchmark_kind=algorithm,
        n=n,
        n_points=n_points,
        backend=meta.backend,
        optimized=meta.optimized,
        cache=meta.cache,
        line_search="safegoldenratio",
        iterations=algorithm == "fw" ? Int(cfg["fw_iterations"]) : Int(cfg["afw_iterations"]),
        elapsed_s=NaN,
        alloc_bytes=0,
        alloc_count=0,
        final_primal=NaN,
        checksum=NaN,
        status="failed",
        error_message=message,
    )
end

function write_output(output_file::AbstractString, row)
    mkpath(dirname(output_file))
    CSV.write(output_file, DataFrame([row]))
    println("Wrote $(output_file)")
end

function main(args::Vector{String})
    cfg_path, algorithm, n, backend_tag = parse_args(args)
    cfg = load_config(cfg_path)
    out_dir = results_dir(cfg_path, cfg)
    output_file = joinpath(out_dir, "csv", "$(case_id(algorithm, n, backend_tag)).csv")
    vertices = generate_point_clouds(cfg, n)

    try
        workload = () -> run_probe(cfg, algorithm, vertices, backend_tag)
        measured, elapsed_s, alloc_bytes, alloc_count = measure(workload)
        result, iterations = measured
        row = success_row(cfg, algorithm, n, backend_tag, elapsed_s, alloc_bytes, alloc_count, result, iterations)
        write_output(output_file, row)
    catch err
        row = failure_row(cfg, algorithm, n, backend_tag, sprint(showerror, err))
        write_output(output_file, row)
        rethrow(err)
    end
end

try
    main(ARGS)
catch err
    showerror(stderr, err)
    println(stderr)
    println(stderr, usage())
    exit(1)
end
