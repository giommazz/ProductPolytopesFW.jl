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
    using GLPK

    # Plotting utilities needed to run FrankWolfe/plot_utils.jl
    include(joinpath(dirname(pathof(FrankWolfe)), "../examples/plot_utils.jl"))
    include("config.jl")
    include("objective_functions.jl")
    include("lmo_utils.jl")
    include("product_algorithms.jl")
    include("utils.jl")
    include("polytopes.jl")

    # From `config.jl`
    export Config, update_config
    # From `objective_functions.jl`
    export objective, gradient!
    # From `lmos.jl`
    export create_product_lmo, find_starting_point, get_lmos
    # From `product_algorithms_algorithms.jl`
    export run_FW, get_solutions
    # From `utils.jl`
    export unique_combinations, generate_rand_float_vector, extract_n_k_iters_from_filename
    # From `polytopes.jl`
    export generate_intersecting_polytopes, save_intersecting_polytopes, load_intersecting_polytopes
    #  from `plot_utils.jl`
    export plot_trajectories

end # module BPCGProduct