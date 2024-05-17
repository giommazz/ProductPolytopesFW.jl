# `utils.jl`
function unique_combinations(config::Config, list::Vector{T}) where T
    # Generate all unique combinations of length k
    return collect(combinations(list, config.k))
end

function generate_rand_float_vector(config::Config; lb=0::T, ub=100::T, seed=42::Int)
    # Set the seed for reproducibility
    Random.seed!(config.seed)
    # Generate random Float64 in [a, b]
    return lb .+ (ub - lb) .* rand(Float64, config.n)
end



TODO: FROM HERE
# Function to extract n and k from the filename
function extract_n_k(filename::String)
    # Regex to match the pattern in the filename
    match = match(r"_n(\d+)_k(\d+)_", filename)
    if match === nothing
        error("Filename format is incorrect.")
    end
    
    # Extract n and k from the matched groups
    n = parse(Int, match.captures[1])
    k = parse(Int, match.captures[2])
    
    return n, k
end