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

function generate_filename(config::Config)
    timestamp = Dates.format(now(), "yyyymmddHHMMSS")
    
    oracle = config.cvxhflag ? "cvxho" : "lmo"
    anc = config.anc_flag ? "anc" : "vert"
    seed = config.seed
    return "k$(config.k)_n$(config.n)_s$(seed)_$(oracle)_$(anc)_t$timestamp"
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
        error("The number of labels ($(length(labels))) does not correspond to the number of FW algorithms run ($(length(traj_data))).\nPlease fix this in your code.")
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

# `pad_log_data` takes `trajectories` (array)
# Input:    `trajectories` contains as many Vectors as the n. of FW variants that have been run
#           each Vector, corresponding to one FW variant run, contains as many 5-entry Tuples as the n. iterations in that run
function pad_log_data(trajectories::Vector{Vector{Tuple{T, T}}}) where T<:Number
    # Find max number of rows among all trajectories
    max_length = maximum(length.(trajectories))
    min_length = minimum(length.(trajectories))
    # For each trajectory, if shorter than `max_length`: append (pad) its last row until length == `max_length`
    padded_data = [length(traj) == max_length ? traj : vcat(traj, fill(traj[end], max_length - length(traj))) for traj in trajectories]
    
    return padded_data, max_length, min_length
end

# Aggregate padded trajectories from multiple FW variants into a DataFrame and write it into a CSV file
function save_padded_logdata_to_csv(
    padded_trajectories::Vector{Vector{NTuple{5, T}}},
    max_length::Int64,
    labels::Vector{String},
    logs_dir::String,
    basename::String
    ) where T <: Number

    # Check that the number of labels matches the number of trajectories
    if length(labels) ≠ length(padded_trajectories)
        error("The number of labels ($(length(labels))) does not correspond to the number of FW algorithms run ($(length(padded_trajectories))).\nPlease fix this in your code.")
    end

    # Transform labels by deleting dashes.
    labels_logs = [replace(l, "-" => "") for l in labels]
    n_variants = length(padded_trajectories)

    # Prepare empty DataFrame. Labels: 'iter' + for each FW variant, 4 columns 'FWVar_pgap', 'FWVar_dual', 'FWVar_dgap', 'FWVar_time'
    df = DataFrame()
    
    # Create or initialize a column named `:iter` in the DataFrame
    # indexing with `!`, we access the column by reference ("inplace" update) → modify actual column rather than a copy
    # using `:` before `iter` defines it as a Symbol, ← column names in Julia DataFrames are typically Symbols
    # `Int64[]` builds the `:iter` column: an empty vector of type Int64
    df[!, :iter] = Int64[]
    
    # For every label in `labels_logs`, create 4 new columns
    for l in labels_logs
        # `Symbol(l * "_pgap")` concatenates `l` with "_pgap" and converts the result into a Symbol, which is then used as column name
        # `Vector{T}()` initializes an empty vector with element type T for that column
        df[!, Symbol(l * "_pgap")] = Vector{T}()
        df[!, Symbol(l * "_dual")] = Vector{T}()
        df[!, Symbol(l * "_dgap")] = Vector{T}()
        df[!, Symbol(l * "_time")] = Vector{T}()
    end

    # Build DataFrame rows, each corresponding to one iteration
    # ------------------------------
    for iter in 1:max_length
        # Initialize Dict to temporarily hold data for current row, keys will be Symbols corresponding to DataFrame column names
        row_dict = Dict{Symbol, Any}()
        row_dict[:iter] = iter
        # For each iteration, extract 4-tuple from each FW variant: ignore 1ˢᵗ element (iter) and retrieve (pgap, dual, dgap, time)
        for fw_variant_i in 1:n_variants
            # Retrieve 5-tuple (iter, pgap, dual, dgap, time)
            tup = padded_trajectories[fw_variant_i][iter]
            fw_variant_baselabel = labels_logs[fw_variant_i]
            # Build keys for each of the four columns
            key_pgap = Symbol(fw_variant_baselabel * "_pgap")
            key_dual = Symbol(fw_variant_baselabel * "_dual")
            key_dgap = Symbol(fw_variant_baselabel * "_dgap")
            key_time = Symbol(fw_variant_baselabel * "_time")
            
            # Assign the corresponding tuple elements to the proper keys in the row dictionary.
            row_dict[key_pgap] = tup[2]
            row_dict[key_dual] = tup[3]
            row_dict[key_dgap] = tup[4]
            row_dict[key_time] = tup[5]
        end
        
        # Push current `row_dict` into DataFrame: add a new row to the DataFrame via `push!`
        # `;` is used to create a NamedTuple from the `row_dict`
        #  `...` is the "splat" operator: unpacks the contents of `row_dict` as keyword arguments into the NamedTuple constructor
        #       → is a concise way to convert a Dict of key=>value pairs into a NamedTuple
        push!(df, (; row_dict...))
    end

    # Write DataFrame to a CSV file.
    CSV.write(joinpath(logs_dir, basename * ".csv"), df)
end