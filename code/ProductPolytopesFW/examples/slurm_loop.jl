#!/usr/bin/env julia
# Generate per-run config/script variants and submit a `(k, n)` grid to Slurm
# through `examples/slurm_experiments.sh`, using hardcoded settings below.
using ProductPolytopesFW

const CONFIG_LINE_TEMPLATE = "config = Config(\"examples/config.yml\")"

# -----------------------------------------------------------------------------
# Hardcoded run settings
# -----------------------------------------------------------------------------
const SCRIPT_TEMPLATE = "examples/compute_intersection_vertex_facet_full_warmup.jl"
const RESULTS_DIR = "examples/results_linesearch_vertex_facet"
const BASE_CONFIG = "examples/config.yml"
const K_VALUES = [2]
const N_VALUES = [10000]
const SEED_START = 42
const DRY_RUN = false

"""
    validate_hardcoded_settings()

Validate the hardcoded run settings and return normalized absolute paths/values.
"""
function validate_hardcoded_settings()
    script_template = abspath(SCRIPT_TEMPLATE)
    results_dir = abspath(RESULTS_DIR)
    base_config = abspath(BASE_CONFIG)

    isfile(script_template) || error("Script template not found: $script_template")
    isfile(base_config) || error("Base config not found: $base_config")
    (endswith(base_config, ".yml") || endswith(base_config, ".yaml")) ||
        error("Base config must be a .yml/.yaml file: $base_config")

    script_kind = if occursin("vertex_facet", lowercase(script_template))
        "vertex_facet"
    elseif occursin("point_clouds", lowercase(script_template))
        "point_clouds"
    else
        error("SCRIPT_TEMPLATE must contain either 'vertex_facet' or 'point_clouds': $script_template")
    end
    occursin(script_kind, lowercase(results_dir)) ||
        error("RESULTS_DIR must contain '$script_kind' to match SCRIPT_TEMPLATE. Found: $results_dir")

    isempty(K_VALUES) && error("K_VALUES must contain at least one integer.")
    isempty(N_VALUES) && error("N_VALUES must contain at least one integer.")
    all(k -> k > 0, K_VALUES) || error("All K_VALUES must be positive integers.")
    all(n -> n > 0, N_VALUES) || error("All N_VALUES must be positive integers.")
    SEED_START >= 0 || error("SEED_START must be non-negative.")

    return script_template, results_dir, base_config, Int.(K_VALUES), Int.(N_VALUES), Int(SEED_START), Bool(DRY_RUN)
end

"""
    write_modified_config_for_run(base_config, config_dir, k, n, seed)

Create a per-run YAML config by overriding `(k, n, seed)` in the base config.
The generated filename includes iteration budgets to avoid collisions across reruns.
"""
function write_modified_config_for_run(
    base_config::Config,
    config_dir::AbstractString,
    k::Int,
    n::Int,
    seed::Int,
)
    modified = modify_config(base_config, k=k, n=n, seed=seed)
    run_tag = "k$(k)_n$(n)_i$(modified.max_iterations)_io$(modified.max_iterations_opt)_s$(seed)"
    config_filename = joinpath(config_dir, "config_$(run_tag).yml")
    return write_config(modified, config_filename), run_tag
end

"""
    write_script_variant(script_template, script_dir, config_filename, run_tag)

Create a per-run Julia script that points to the generated config file.
The template is required to contain exactly one `config = Config("...")` line.
"""
function write_script_variant(
    script_template::AbstractString,
    script_dir::AbstractString,
    config_filename::AbstractString,
    run_tag::AbstractString,
)
    src_text = read(script_template, String)
    # Replace exactly one config-loading line so the generated script is tied to one run config.
    occurrences = length(split(src_text, CONFIG_LINE_TEMPLATE)) - 1
    occurrences == 1 || error(
        "Expected exactly one occurrence of '$CONFIG_LINE_TEMPLATE' in $(script_template), found $(occurrences).",
    )

    replacement = "config = Config(\"$(abspath(config_filename))\")"
    new_text = replace(src_text, CONFIG_LINE_TEMPLATE => replacement)

    script_stem = splitext(basename(script_template))[1]
    script_filename = joinpath(script_dir, "$(script_stem)_$(run_tag).jl")
    write(script_filename, new_text)
    return abspath(script_filename)
end

"""
    safe_submit(cmd; dry_run=false)

Submit one `sbatch` command, or print it without submitting when `dry_run=true`.
"""
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

"""
    main()

Generate per-run config/script files for each hardcoded `(k, n)` pair and submit
the corresponding Slurm jobs.
"""
function main()
    script_template, results_dir, base_config_filename, k_values, n_values, seed_start, dry_run = validate_hardcoded_settings()
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
        run_config, run_tag = write_modified_config_for_run(base_config, generated_config_dir, k, n, seed)
        run_script = write_script_variant(script_template, generated_script_dir, run_config, run_tag)

        # Keep the Slurm CPU request synchronized with the current `k` value.
        slurm_cmd = `sbatch --cpus-per-task=$(k) $(slurm_wrapper) $(run_script) $(results_dir) $(run_config)`
        println("Submitting job k=$(k), n=$(n), seed=$(seed)")
        safe_submit(slurm_cmd; dry_run=dry_run)
        seed += 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end