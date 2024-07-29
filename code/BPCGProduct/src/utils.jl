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

function generate_filename(config::Config) where T
    timestamp = Dates.format(now(), "yyyymmddHHMMSS")
    
    oracle = config.cvxhflag ? "cvxho" : "lmo"
    anc = config.anc_flag ? "anc" : "vert"
    return "n$(config.n)_k$(config.k)_$(oracle)_$(anc)_t$timestamp"
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

# Compute avg iteration time and tot time of FW algo execution
# Call as     `log_data(trajectories_i, ["C-BC-FW", "F-BC-BPCG"], "filename")`
function log_data(traj_data::Vector{Any}, labels::Vector{String}, basename::String)
    
    if length(labels) ≠ length(traj_data)
        error("There are $(length(labels)) labels but only $(length(traj_data)) FW algorithms were executed")
    end

    # Store results in DataFrame object
    times = DataFrames.DataFrame(label=String[], iters=Int64[], tot_time=Float64[], avg_time=Float64[])

    for i in 1:length(traj_data)
        label = labels[i]
        # Compute total time summing time over all iterations
        tot_time = sum(iter -> iter[5], traj_data[i])
        iters = length(traj_data[i])
        avg_time = tot_time / iters
        
        println("Label: $label, Total Time: $tot_time, Average Time per Iteration: $avg_time")
        push!(times, (label, iters, tot_time, avg_time))    
    end
    CSV.write(basename*".csv", times)
end