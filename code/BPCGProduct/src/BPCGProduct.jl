module BPCGProduct
    using FrankWolfe
    using Plots
    using Random
    using Combinatorics
    using YAML


    include(joinpath(dirname(pathof(FrankWolfe)), "../examples/plot_utils.jl"))

    include("objective_functions.jl")
    include("lmos.jl")
    include("fw_algorithms.jl")
    include("utils.jl")
    include("config.jl")

    # From `objective_functions.jl`
    export f, grad!
    # From `lmos.jl`
    export create_product_lmo, find_starting_point
    # From `fw_algorithms.jl`
    export run_FW
    # From `utils.jl`
    export unique_combinations, generate_rand_float_vector
    # From `config.jl`
    export Config


end # module BPCGProduct