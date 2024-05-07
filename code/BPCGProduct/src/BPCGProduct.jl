# `BPCGProduct.jl`
module BPCGProduct
    using FrankWolfe
    using Plots
    using Random
    using Combinatorics
    using YAML
    using JuMP
    using LinearAlgebra
    using Plots
    using MathOptInterface
    using CDDLib
    using Polyhedra
    using Makie

    include(joinpath(dirname(pathof(FrankWolfe)), "../examples/plot_utils.jl"))

    include("objective_functions.jl")
    include("lmos.jl")
    include("fw_algorithms.jl")
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
    export generate_polytope, generate_pyramid_polytope, plot_polytope_2d
    # From `Polyhedra`
    export polyhedron

end # module BPCGProduct