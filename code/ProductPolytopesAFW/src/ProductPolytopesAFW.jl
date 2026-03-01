# `ProductPolytopesAFW.jl`
module ProductPolytopesAFW
    using FrankWolfe
    using Plots           # Plotting utilities needed to run FrankWolfe/plot_utils.jl
    using Random
    using Combinatorics
    using YAML
    using JuMP
    using LinearAlgebra
    using CDDLib
    using Polyhedra
    using FileIO
    using Dates
    using GLPK, SCIP, HiGHS
    using CSV, DataFrames
    using LinearAlgebra          # brings dot and BLAS.axpy!
    using Statistics

    # Plotting utilities needed to run FrankWolfe/plot_utils.jl
    include(joinpath(dirname(pathof(FrankWolfe)), "../examples/plot_utils.jl"))
    include("config.jl")
    include("frankwolfe_patches.jl")
    include("objective_functions.jl")
    include("lmo_utils.jl")
    include("utils.jl")
    include("polytope_utils.jl")
    include("product_algorithms.jl")
    include("polytope_generation.jl")
    include("plotting_utils.jl")

    # From `config.jl`
    export Config, print_config, modify_config, write_config
    # From `objective_functions.jl`
    export convex_feasibility_objective, convex_feasibility_gradient!
    # From `lmos.jl`
    export create_product_lmo, find_starting_point, create_lmos
    # From `product_algorithms_algorithms.jl`
    export run_BlockCoordinateFW, run_FullFW, run_FullAFW, run_AlternatingProjections, AwayStep
    # From `frankwolfe_patches.jl`
    export set_fw_weight_purge_default_override!, fw_weight_purge_default_override
    # From `utils.jl`
    export ensure_dir, unique_combinations, generate_rand_float_vector, extract_n_k_from_filename, base_name, 
        approxequal, log_times, pad_log_data, save_logdata_to_csv,
        push_to_trajectories!, save_trajectories, load_trajectories, compute_primal_gap, best_seen_solution, get_k_n_from_logstring
    # From `polytope_generation.jl` / `polytope_utils.jl`
    export generate_polytopes, compute_distance, save_polytopes, load_polytopes, generate_filename
    # From `plot_utils.jl`
    export plot_trajectories
    # From `plotting_utils.jl`
    export plot_time_only, cutoff_log_shortest_time, cutoff_time, load_fw_trajectories_i, load_fw_trajectories_ni, avg_over_logs

end # module ProductPolytopesAFW
