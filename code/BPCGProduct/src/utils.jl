# `utils.jl`
function unique_combinations(list, config)
    # Generate all unique combinations of length k
    return collect(combinations(list, config.k))
end

function generate_rand_float_vector(config, lb=0, ub=100, seed=42)
    # Set the seed for reproducibility
    Random.seed!(config.seed)
    # Generate random Float64 in [a, b]
    return lb .+ (ub - lb) .* rand(Float64, config.n)
end