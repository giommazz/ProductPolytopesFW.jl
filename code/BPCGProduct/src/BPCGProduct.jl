# `BPCGProduct.jl`
module BPCGProduct
    using FrankWolfe
    # using Plots           # Plotting utilities needed to run FrankWolfe/plot_utils.jl
    using Random
    using Combinatorics
    using YAML
    using JuMP
    using LinearAlgebra
    using MathOptInterface
    using CDDLib
    using Polyhedra
    using Makie
    using GLMakie
    using GLPK

    # Plotting utilities needed to run FrankWolfe/plot_utils.jl
    # include(joinpath(dirname(pathof(FrankWolfe)), "../examples/plot_utils.jl"))

    include("objective_functions.jl")
    include("lmo_utils.jl")
    include("product_algorithms.jl")
    include("utils.jl")
    include("config.jl")
    include("polytopes.jl")

    # From `objective_functions.jl`
    export objective, gradient!
    # From `lmos.jl`
    export create_product_lmo, find_starting_point
    # From `fw_algorithms.jl`
    export run_FW
    # From `utils.jl`
    export unique_combinations, generate_rand_float_vector
    # From `config.jl`
    export Config
    # From `polytopes.jl`
    export generate_polytope, find_vertex_in_polytope, generate_simplex, setup_translated_polytope_B

end # module BPCGProduct