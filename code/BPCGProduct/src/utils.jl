# `utils.jl`
function unique_combinations(config::Config, list::Vector)
    # Generate all unique combinations of length k
    return collect(combinations(list, config.k))
end

function generate_rand_float_vector(config::Config, lb=0::T, ub=100::T, seed=42::Int)
    # Set the seed for reproducibility
    Random.seed!(config.seed)
    # Generate random Float64 in [a, b]
    return lb .+ (ub - lb) .* rand(Float64, config.n)
end