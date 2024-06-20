# `utils.jl`

# Generate all unique combinations of length k and convert to a list
function unique_combinations(config::Config, list::Vector{T}) where T
    return collect(combinations(list, config.k))
end

function generate_rand_float_vector(config::Config; lb=0, ub=100, seed=42::Int)
    # Set the seed for reproducibility
    Random.seed!(config.seed)
    # Generate random Float64 in [a, b]
    return lb .+ (ub - lb) .* rand(Float64, config.n)
end

function generate_filename(config::Config, vertices::Vector{Matrix{T}}) where T
    sizes = [size(poly_vertices)[1] for poly_vertices in vertices]
    timestamp = Dates.format(now(), "yyyymmddHHMMSS")
    return "intersecting_polytopes_n$(config.n)_k$(config.k)_v$(join(sizes, "-"))_t$timestamp.jld2"
    return "intersecting_polytopes_n$(config.n)_k$(config.k)_mi$(config.max_iterations)_v$(join(sizes, "-"))_t$timestamp.jld2"

end

# Function to extract `n`, `k`, and `max_iterations` from the filename
function extract_n_k_from_filename(filename::String)
    # Regex to match the pattern in the filename
    m = match(r"_n(\d+)_k(\d+)_", filename)
    if m === nothing
        error("Filename format is incorrect.")
    end
    
    # Extract `n`, `k`, and `max_iterations` from the matched groups
    n = parse(Int, m.captures[1])
    k = parse(Int, m.captures[2])
    
    return n, k
end

# Takes any variable and returns its name as a String. Use as a macro: @var_name
macro var_name(var)
    return string(var)
end

# Extract base name without extension, for files of following type: "intersecting_polytopes_n5_k2_v10-11_t20240524121720.jld2"
function base_name(filename::String)
    return splitext(filename)[1]
end

# Check if two vectors/numbers are equal (i.e., their difference is equal at each element), up to given tolerance
# Return true if yes 
function approxequal(a, b; tol=1e-04) return all(abs.(a .- b) .≤ tol) end