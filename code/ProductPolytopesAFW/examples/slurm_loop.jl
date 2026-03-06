#!/usr/bin/env julia
using ProductPolytopesAFW

const CONFIG_LINE_TEMPLATE = "config = Config(\"examples/config.yml\")"

function usage()
    return """
Usage:
  julia --project=. examples/slurm_loop.jl <script_template> <results_dir> <base_config> <k_csv> <n_csv> <seed_start> [--dry-run]

Arguments:
  <script_template>  Experiment script template (point-cloud or vertex-facet script).
  <results_dir>      Output directory passed to slurm_experiments.sh.
  <base_config>      Base YAML config file.
  <k_csv>            Comma-separated k values (example: 2,3,4).
  <n_csv>            Comma-separated n values (example: 103,207).
  <seed_start>       Initial integer seed, incremented after each submitted job.
  --dry-run          Optional flag: print sbatch commands without submitting.
"""
end

function parse_int_csv(csv::AbstractString, field_name::AbstractString)
    values = Int[]
    for token in split(csv, ",")
        token = strip(token)
        isempty(token) && error("Invalid $field_name: empty value in '$csv'.")
        value = try
            parse(Int, token)
        catch
            error("Invalid $field_name: '$token' is not an integer.")
        end
        value > 0 || error("Invalid $field_name: all values must be positive integers.")
        push!(values, value)
    end
    isempty(values) && error("Invalid $field_name: no values parsed.")
    return values
end

function parse_seed(seed_raw::AbstractString)
    seed = try
        parse(Int, seed_raw)
    catch
        error("Invalid seed_start '$seed_raw': must be an integer.")
    end
    seed >= 0 || error("Invalid seed_start '$seed_raw': must be non-negative.")
    return seed
end

function parse_args(args::Vector{String})
    if !(length(args) in (6, 7))
        error(usage())
    end
    dry_run = false
    if length(args) == 7
        args[7] == "--dry-run" || error("Unknown option '$(args[7])'.\n\n$(usage())")
        dry_run = true
    end

    script_template = abspath(args[1])
    results_dir = abspath(args[2])
    base_config = abspath(args[3])
    k_values = parse_int_csv(args[4], "k_csv")
    n_values = parse_int_csv(args[5], "n_csv")
    seed_start = parse_seed(args[6])

    isfile(script_template) || error("Script template not found: $script_template")
    isfile(base_config) || error("Base config not found: $base_config")
    (endswith(base_config, ".yml") || endswith(base_config, ".yaml")) || error("Base config must be a .yml/.yaml file: $base_config")

    return script_template, results_dir, base_config, k_values, n_values, seed_start, dry_run
end

function write_modified_config_for_run(
    base_config::Config,
    config_dir::AbstractString,
    k::Int,
    n::Int,
    seed::Int,
)
    modified = modify_config(base_config, k=k, n=n, seed=seed)
    config_filename = joinpath(config_dir, "config_k$(k)_n$(n)_s$(seed).yml")
    return write_config(modified, config_filename)
end

function write_script_variant(
    script_template::AbstractString,
    script_dir::AbstractString,
    config_filename::AbstractString,
    k::Int,
    n::Int,
    seed::Int,
)
    src_text = read(script_template, String)
    occurrences = length(split(src_text, CONFIG_LINE_TEMPLATE)) - 1
    occurrences == 1 || error(
        "Expected exactly one occurrence of '$CONFIG_LINE_TEMPLATE' in $(script_template), found $(occurrences).",
    )

    replacement = "config = Config(\"$(abspath(config_filename))\")"
    new_text = replace(src_text, CONFIG_LINE_TEMPLATE => replacement)

    script_stem = splitext(basename(script_template))[1]
    script_filename = joinpath(script_dir, "$(script_stem)_k$(k)_n$(n)_s$(seed).jl")
    write(script_filename, new_text)
    return abspath(script_filename)
end

function safe_submit(cmd::Cmd; dry_run::Bool=false)
    if dry_run
        println("[dry-run] ", cmd)
        return
    end
    try
        run(cmd)
    catch err
        @warn "slurm submission failed" cmd=cmd error=err
    end
end

function main(args::Vector{String})
    script_template, results_dir, base_config_filename, k_values, n_values, seed_start, dry_run = parse_args(args)
    base_config = Config(base_config_filename)

    generated_root = joinpath(results_dir, "slurm_generated")
    generated_config_dir = joinpath(generated_root, "configs")
    generated_script_dir = joinpath(generated_root, "scripts")
    mkpath(generated_config_dir)
    mkpath(generated_script_dir)
    mkpath(results_dir)

    slurm_wrapper = abspath(joinpath(dirname(@__FILE__), "slurm_experiments.sh"))
    isfile(slurm_wrapper) || error("SLURM wrapper script not found: $slurm_wrapper")

    seed = seed_start
    for k in k_values, n in n_values
        run_config = write_modified_config_for_run(base_config, generated_config_dir, k, n, seed)
        run_script = write_script_variant(script_template, generated_script_dir, run_config, k, n, seed)

        slurm_cmd = `sbatch --cpus-per-task=$(k) $(slurm_wrapper) $(run_script) $(results_dir) $(run_config)`
        println("Submitting job k=$(k), n=$(n), seed=$(seed)")
        safe_submit(slurm_cmd; dry_run=dry_run)
        seed += 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main(ARGS)
    catch err
        println(stderr, err)
        println(stderr, "")
        println(stderr, usage())
        exit(1)
    end
end
