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
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    
    oracle = config.cvxhflag ? "cvxho" : "lmo"
    anc = config.anc_flag ? "anc" : "vert"
    seed = config.seed
    return "k$(config.k)_n$(config.n)_i$(config.max_iterations)_s$(seed)_$(oracle)_$(anc)_t$timestamp"
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
function log_times(traj_data::Vector{Any}, labels::Vector{String}, basename::String)
    
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
function pad_log_data(trajectories::Vector{Any})
    # Find max number of rows among all trajectories
    max_length = maximum(length.(trajectories))
    min_length = minimum(length.(trajectories))
    # For each trajectory, if shorter than `max_length`: append (pad) its last row until length == `max_length`
    padded_data = [length(traj) == max_length ? traj : vcat(traj, fill(traj[end], max_length - length(traj))) for traj in trajectories]
    
    return padded_data, max_length, min_length
end

# Aggregate padded trajectories from multiple FW variants into a DataFrame and write it into a CSV file
function save_logdata_to_csv(
    padded_trajectories::Vector{Vector{Any}},
    max_length::Int64,
    labels::Vector{String},
    logs_dir::String,
    basename::String
    )

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
        df[!, Symbol(l * "_pgap")] = Vector{Float64}()
        df[!, Symbol(l * "_dual")] = Vector{Float64}()
        df[!, Symbol(l * "_dgap")] = Vector{Float64}()
        df[!, Symbol(l * "_time")] = Vector{Float64}()
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


"""
    load_fw_trajectories(path::String; T=Float64)

Input: .csv log at `path` of the form "iter,FW1_pgap,FW1_dual,FW1_dgap,FW1_time,FW2_pgap, ..."
Output: Vector{Vector{NTuple{5, T}}} of trajectories such that each NTuple is (iter, pgap, dual, dgap, time).
"""
function load_fw_trajectories(path::String)
    
    # `splitext` returns (root, extension)
    ext = splitext(path)[2]
    if lowercase(ext) != ".csv"
        throw(ArgumentError("`path` must end in \".csv\" but extension is \"$ext\""))
    end

    # read .csv file
    df = CSV.read(path, DataFrame)
    
    # determine n. of FW variants that were executed (first column is :iter)
    ncols = ncol(df)
    @assert ncols ≥ 5 "Expected ≥5 columns (iter + at least one variant's 4 columns)"
    n_variants = Int((ncols - 1) ÷ 4)

    # pre-allocate the output: one vector per variant
    trajectories = Vector{Vector{Any}}(undef, n_variants)
    for i in 1:n_variants
        trajectories[i] = Vector{Any}(undef, nrow(df))
    end

    # extract column names
    cols = names(df)

    # define allowed tokens
    FW_TOKENS = ["AP", "BPFW", "AFW", "FW", "BC", "C", "F"]
    # sort tokens by descending length so that multi‑letter tokens (e.g. "BPFW") match before shorter ones ("FW")
    sorted_tokens = sort(FW_TOKENS, by=length, rev=true)
    # join the sorted tokens with "|" to build the alternation part of a regex pattern
    token_pattern = join(sorted_tokens, "|")
    # compile the final regex that will match any one of the tokens
    token_re = Regex("($token_pattern)")

    # prepare empty array of strings to hold the final human-readable labels
    variant_labels = String[]
    # for each FW variant (each group of 4 columns in the .csv logfile)...
    for v in 1:n_variants
        # a) pick header name for this variant's _pgap column,b ) convert to plain string, c) split on "_" to isolate block name (e.g. "CBCFW")
        base_name = split(String(cols[1 + (v-1)*4 + 1]), "_")[1]  # e.g. "CBCFW"
        # find all occurrences of tokens in that block name, yielding ["C","BC","FW"]
        parts = [m.match for m in eachmatch(token_re, base_name)]
        # join the parts with "-" to get "C-BC-FW", then append to the labels array
        push!(variant_labels, join(parts, "-")) # e.g. "C-BC-FW"
    end

    # row indices
    for row_idx in 1:nrow(df)
        # get only columns ≠ :iter
        iter = Int64(df[row_idx, cols[1]])
        # for each variant, grab its 4 fields
        for FWvar_index in 1:n_variants
            # column indices: compute "starting column -1" for each each variant...
            base = 1 + (FWvar_index - 1)*4
            # ...then retrieve tuple indices for each variant and each row
            pgap_col = cols[base+1]
            dual_col = cols[base+2]
            dgap_col = cols[base+3]
            time_col = cols[base+4]

            # now you have row indices and column indices, convert values from String to Number
            pgap = Float64(df[row_idx, pgap_col])
            dual = Float64(df[row_idx, dual_col])
            dgap = Float64(df[row_idx, dgap_col])
            time = Float64(df[row_idx, time_col])

            trajectories[FWvar_index][row_idx] = (iter, pgap, dual, dgap, time)
        end
    end

    return trajectories, variant_labels
end


# Push trajectory data into appropriate vector (intersecting or non-intersecting) based on the flag `ni_flag` 
function push_to_trajectories!(
    ni_flag::Bool,
    trajectory_data_curr::Vector{Any},
    trajectories_ni::Vector{Any}, # initially empty vector, will be filled with the ni_trajectory data
    trajectories_i::Vector{Any}, # initially empty vector, will be filled with the i_trajectory data
    primal::Float64
    )
    # `primal` ≠ 0.0: the polytopes don't intersect
    if ni_flag
        # Replace "Primal" with "Primal Gap" in the FW log, i.e., replace f(x) with f(x) - `primal` 
        trajectory_data_pg = compute_primal_gap(trajectory_data_curr, primal)
        push!(trajectories_ni, trajectory_data_pg)

    # `primal` == 0.0: the polytopes do intersect
    else    
        push!(trajectories_i, trajectory_data_curr)
    end
end
# (Multiple dispatch)
function push_to_trajectories!(
    ni_flag::Bool,
    trajectory_data_curr::Vector{Any},
    trajectories::Vector{Any},
    primal::Float64
    )

    # `primal` ≠ 0.0: the polytopes don't intersect
    if ni_flag
        # Replace "Primal" with "Primal Gap" in the FW log, i.e., replace f(x) with f(x) - `primal` 
        trajectory_data_pg = compute_primal_gap(trajectory_data_curr, primal)
        push!(trajectories, trajectory_data_pg)
    # `primal` == 0.0: the polytopes do intersect
    else    
        push!(trajectories, trajectory_data_curr)
    end
end

# Save trajectory data to given .jld2 file
function save_trajectories(
    filename::String,
    trajectories::Vector{Any}
    )
    
    save(filename, Dict("trajectories" => trajectories))
    println("Saving data to $filename")
end
# (Multiple dispatch) handle several trajectory data elements
function save_trajectories(
    filename::String,
    trajectories...
    )
    
    dict = Dict{String, Vector{Any}}()
    for td in trajectories
        dict[@var_name(td)] = td
    end
    save(filename, dict)
    println("Saving data to $filename with multiple trajectory data")
end

# Load data from .jld2 file
function load_trajectories(filename::String)
    
    # Load the file
    f = load(filename) 
    # Extract the keys and values from the loaded dictionary
    trajectory_keys = keys(f)
    trajectory_values = [f[k] for k in trajectory_keys]
    
    println("Loaded data from $filename with keys: $keys")
    
    return trajectory_values
end

function compute_primal_gap(trajectories_curr::Vector{Any}, opt::Float64)
    trajectories_curr_pg = deepcopy(trajectories_curr)
    # Iterate over each tuple within the current block
    for i in 1:length(trajectories_curr_pg)
        # for FW algorithms:            iter, primal, dual,           dgap, time
        # for Alternating projections:  iter, infeas, partial infeas, dgap, time
        iter, primal, dual, dgap, time = trajectories_curr_pg[i]    
        # Compute primal gap
        pgap = primal - opt
        # Replace current tuple with a new one including the primal gap instead of the primal value
        trajectories_curr_pg[i] = (iter, pgap, dual, dgap, time)
    end
    neg_idx = findall(<(0.0), [t[2] for t in trajectories_curr])
    lni = length(neg_idx)
    if lni > 0
        error("There are $(lni) negative values when computing the primal gap. Check `plotting_utils.jl`")
    end
    return trajectories_curr_pg
end