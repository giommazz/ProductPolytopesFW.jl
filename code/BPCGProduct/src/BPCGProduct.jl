# `BPCGProduct.jl`
module BPCGProduct
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

    # Plotting utilities needed to run FrankWolfe/plot_utils.jl
    include(joinpath(dirname(pathof(FrankWolfe)), "../examples/plot_utils.jl"))
    include("config.jl")
    include("objective_functions.jl")
    include("lmo_utils.jl")
    include("product_algorithms.jl")
    include("utils.jl")
    include("polytopes.jl")
    include("plotting_utils.jl")

    # From `config.jl`
    export Config, print_config
    # From `objective_functions.jl`
    export objective, gradient!
    # From `lmos.jl`
    export create_product_lmo, find_starting_point, create_lmos
    # From `product_algorithms_algorithms.jl`
    export run_FW
    # From `utils.jl`
    export unique_combinations, generate_rand_float_vector, extract_n_k_from_filename
    # From `polytopes.jl`
    export generate_polytopes, save_polytopes, load_polytopes
    # From `plot_utils.jl`
    export plot_trajectories
    # From `plotting_utils.jl`
    export print_trajdata, compute_primal_gap

end # module BPCGProduct