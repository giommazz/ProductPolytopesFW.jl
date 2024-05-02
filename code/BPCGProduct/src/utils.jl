function unique_combinations(list, k)
    # Generate all unique combinations of length k
    combos = combinations(list, k)
    # Include repeated element combinations
    repeated_combos = [fill(entry, k) for entry in list]
    # Concatenate and return unique combinations
    unique_combos = vcat(collect(combos), repeated_combos)
    return unique_combos
end

function generate_rand_float_vector(lb=0, ub=100, seed=42)
    # Set the seed for reproducibility
    Random.seed!(seed)
    # Generate random Float64 in [a, b]
    return lb .+ (ub - lb) .* rand(Float64, n)
end