# `plotting_utils.jl`

"""
Given a vector of trajectories (each being a vector of 5‑tuples), extract time-based info and plot two subplots:
    - left: pgap (tuple field 2) over time (tuple field 5)
    - right: dgap (FW gap) (tuple field 4) over time

Arguments:
- `trajectories`: contains data to be plotted
- `labels`: as many as `length(trajectories)`
- `xscalelog`, `yscalelog`: axis scaling
- `primal_offset`: vertical offset for primal subplot (if needed).
- `offset`: starting index in each trajectory (default 2, as in the FrankWolfe.jl package)

Return: Plot object with two subplots
"""
function plot_time_only(
    config::Config,
    trajectories::Vector{Vector{Any}},
    labels::Vector{String};
    xscalelog::Bool = true,   # Default to true for log scale on x-axis
    yscalelog::Bool = true,   # Default to true for log scale on y-axis
    legend_position = :bottomleft,
    lstyle = fill(:solid, length(trajectories)),
    line_width = 1.3,
    primal_offset = 1e-8,
    offset = 2,
    size::Tuple{Int, Int} = (1200, 400)
    )
    # Custom colorblind palette from https://www.color-hex.com/color-palette/49436
    colorblind_palette = [
        "#E69F00",  # Orange
        "#CC79A7",  # Pink
        "#0072B2",  # Blue
        "#009E73",  # Teal
        "#F0E442"   # Yellow
        ]

    legend_lbl = "k=$(config.k), n=$(config.n)"   # text that should appear once
    
    # Create empty plot for the "primal gap" (pgap) over time.
    # We let Plots choose tick positions automatically (instead of hardcoding them)
    plt_primal = plot(
        xaxis = xscalelog ? :log10 : :identity,
        yaxis = yscalelog ? :log10 : :identity,
        xlabel = "Time (s)",
        ylabel = "Primal",
        legend = legend_position,
        legendtitle = legend_lbl,
        lw = line_width,
        size = size,
        xguidefontsize = 12,
        yguidefontsize = 12,
        left_margin = 10Plots.mm        # extra space to visualize the "Primal" label
        )

    # Create empty plot for the "dual gap" (FW gap) over time.
    plt_gap = plot(
        xaxis = xscalelog ? :log10 : :identity,
        yaxis = yscalelog ? :log10 : :identity,
        xlabel = "Time (s)",
        ylabel = "FW gap",
        legend = legend_position,
        legendtitle = legend_lbl,
        lw = line_width,
        size = size
        )
        
    # Loop through each trajectory and add its data.
    # We check for sufficient length so that we don't pass an empty vector to the plotting routines.
    for (i, trajectory) in enumerate(trajectories)
        if length(trajectory) < offset
            continue  # Skip this trajectory if it does not contain enough points.
        end

        # Extract the data from the tuple assuming:
        # - Field 2 is the primal gap (pgap)
        # - Field 4 is the FW gap (dual gap)
        # - Field 5 is time.
        times = [trajectory[j][5] for j in offset:length(trajectory)]
        if isempty(times)
            continue  # Skip if no valid time points have been extracted.
        end
        primal_vals = [trajectory[j][2] + primal_offset for j in offset:length(trajectory)]
        gap_vals = [trajectory[j][4] for j in offset:length(trajectory)]
        
        # Choose a color (cycle through the palette if more series than colors are present)
        color_choice = (length(trajectories) <= length(colorblind_palette)) ?
                          colorblind_palette[i] :
                          colorblind_palette[mod1(i, length(colorblind_palette))]
        
        # Plot the current series on both subplots.
        plot!(plt_primal, times, primal_vals, label = labels[i], linestyle = lstyle[i], lw = line_width, color = color_choice)
        plot!(plt_gap, times, gap_vals, label = labels[i], linestyle = lstyle[i], lw = line_width, color = color_choice)
    end

    final_plot = plot(
        plt_primal, plt_gap;
        layout        = (1, 2),
        size          = (1200, 400),
        left_margin   = 10Plots.mm,
        bottom_margin = 15Plots.mm,
        legend        = legend_position,       # one legend for both
    )

    return final_plot
end

"""
    cutoff_time(trajectories::Vector{Vector{Any}})
Input: `trajectories`, where each vector is a vector of 5-tuples: (iter, pgap, dual, dgap, time)
Output: 
    - min of the final times among all runs
"""
function cutoff_time(trajectories::Vector{Vector{Any}})
    # Decide `cutoff_time`: ∀ Vectors in `trajectories`, extract 5th element (time) of the last tuple, then compute min among all these times
    return minimum([last(traj)[5] for traj in trajectories])
end
"""
    [Multiple dispatch] cutoff_time(logfiles::Vector{String})
Input:
    - `logfiles`: Vector of logfile names
Output:
    - `global_cutoff_time`: min cutoff time over all trajectories in all logfiles
"""
function cutoff_time(logfiles::Vector{String}, wanted_fw_variants::Vector{String})
    
    # for each logfile, compute that log's cutoff time
    global_cutoff_time = Inf
    for file in logfiles
        # load trajectories
        traj_ni, _, _ = load_fw_trajectories_ni(file; wanted_fw_variants=wanted_fw_variants)
        # find cutoff time for the current `traj_ni`
        cutoff_time = cutoff_time(traj_ni)
        # update global 
        global_cutoff_time = cutoff_time < global_cutoff_time ? cutoff_time : global_cutoff_time
    end
    
    return global_cutoff_time
end





"""
    cutoff_log_shortest_time(trajectories::Vector{Vector{Any}})
Input: `trajectories`, where each vector is a vector of 5-tuples: (iter, pgap, dual, dgap, time)
Output: 
    - `cutoff_time`: min of the final times among all runs
    - `cutoff_trajectories`: truncated `trajectories`, s.t. every one only contains entries where `time` ≤ `cutoff_time`
"""
function cutoff_log_shortest_time(trajectories::Vector{Vector{Any}})
    
    # Decide `cutoff_time`: ∀ Vectors in `trajectories`, extract 5th element (time) of the last tuple, then compute min among all these times
    cutoff_time = cutoff_time(trajectories)
    # Create `cutoff_trajectories`, truncated to earliest finish point: ∀ Vectors in `trajectories`, cut out tuples where `time` ≥ `cutoff_time`
    cutoff_trajectories = cutoff_log_shortest_time(trajectories, cutoff_time)
    
    return cutoff_trajectories, cutoff_time
end
"""
    [Multiple Dispatch] cutoff_log_shortest_time(trajectories::Vector{Vector{Any}}, cutoff_time::Float64) 
Input:
    - `trajectories`, where each vector is a vector of 5-tuples: (iter, pgap, dual, dgap, time)
    - `cutoff_time`: min of the final times among all runs
Output: 
    - truncated `trajectories`, s.t. every one only contains entries where `time` ≤ `cutoff_time`
"""
function cutoff_log_shortest_time(trajectories::Vector{Vector{Any}}, cutoff_time::Float64)
    
    # Create `cutoff_trajectories`, truncated to earliest finish point: ∀ Vectors in `trajectories`, cut out tuples where `time` ≥ `cutoff_time`
    return [filter(tuple -> tuple[5] ≤ cutoff_time, traj) for traj in trajectories]
end





"""
    load_fw_trajectories(path::String; wanted_fw_variants::Vector{String}=String[])

Input:
    - .csv log at `path` of the form "iter,FW1_prim,FW1_dual,FW1_dgap,FW1_time,FW2_prim, ..."
    - `wanted_fw_variants` (optional): a list of wanted labels, in case we only want to plot some of the logged FW variants
Output:
    - `trajectories`: Vector{Vector{Any}} of trajectories, each with entries (iter, prim, dual, dgap, time)
    - `retrieved_fw_labels`
    - `opt`: optimal value for that instance
"""
function load_fw_trajectories_ni(path::String; wanted_fw_variants::Vector{String}=String[])
    
    # `splitext` returns (root, extension)
    ext = lowercase(splitext(path)[2])
    @assert ext == ".csv" "Expected a .csv file"
    @assert lowercase([1:3]) == "ni_" "Expected a file pertaining to an non-intersecting instance, instead found $(split(path, "/"))"
    
    # read .csv file
    df = CSV.read(path, DataFrame)
    # extract column names
    cols = names(df)
    # determine n. of FW variants that were executed (columns are :iter,:opt, then 4×n_variants columns)
    ncols = ncol(df)
    @assert ncols ≥ 6 "Expected ≥6 columns (iter, opt, plus at least one variant's 4 columns)"
    n_variants = Int((ncols - 2) ÷ 4)

    
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
        # a) pick header name for this variant's _prim column
        # b) convert to plain string
        # c) split on "_" to isolate block name (e.g. "CBCFW")
        # after iter,opt the first "prim" of current FW variant sits at index = 2 + (v-1)*4 + 1
        prim_header = String(cols[2 + (v-1)*4 + 1])         # e.g. "CBCFW_prim"
        block_name  = split(prim_header, "_")[1]    # e.g. "CBCFW"

        # d) find all occurrences of tokens in that block name, yielding ["C","BC","FW"]
        parts       = [m.match for m in eachmatch(token_re, block_name)]
        # e) join the parts with "-" to get "C-BC-FW", then append to the labels array
        push!(variant_labels, join(parts, "-"))     # yields "C-BC-FW", etc.
    end

    # If user asked for a subset, filter indices + labels
    if isempty(wanted_fw_variants)
        wanted_idxs = 1:n_variants
    else
        # retrieve wanted FW labels from `variant_labels` created above
        wanted_idxs = findall(lbl -> lbl in wanted_fw_variants, variant_labels)
        # check that `idxs` is not empty 
        @assert !isempty(wanted_idxs) "No matching FW variants: $(wanted_fw_variants)"
    end
    retrieved_fw_labels   = variant_labels[wanted_idxs]

    # *****************************
    # extract the (constant) opt value from the second column of row 1
    opt = Float64(df[1, cols[2]])

    # pre-allocate the output: one vector per variant
    trajectories = Vector{Vector{Any}}(undef, length(retrieved_fw_labels))
    for i in 1:length(retrieved_fw_labels)
        trajectories[i] = Vector{Any}(undef, nrow(df))
    end

    # row indices
    for row_idx in 1:nrow(df)
        
        # read iteration number
        iter = Int64(df[row_idx, cols[1]])
        # for each selected FW variant, grab its 4 fields
        for local_idx in 1:length(wanted_idxs)
            orig_idx = wanted_idxs[local_idx]
            # skip the first two columns (iter,opt), then block of 4 per FW variant
            base     = 2 + (orig_idx - 1) * 4
            # now you have `row_idx` and column indices, convert values from String to Number
            prim = Float64(df[row_idx, cols[base + 1]])
            dual = Float64(df[row_idx, cols[base + 2]])
            dgap = Float64(df[row_idx, cols[base + 3]])
            time = Float64(df[row_idx, cols[base + 4]])

            trajectories[local_idx][row_idx] = (iter, prim, dual, dgap, time)
        end
    end

    return trajectories, retrieved_fw_labels, opt
end





function load_fw_trajectories_i(path::String; wanted_fw_variants::Vector{String}=String[])
    # ---------------- basic checks -------------------------------------------------
    @assert lowercase(splitext(path)[2]) == ".csv" "Expected a .csv file"
    @assert lowercase([1:2]) == "i_" "Expected a file pertaining to an intersecting instance, instead found $(split(path, "/"))"

    df   = CSV.read(path, DataFrame)
    cols = names(df)

    ncols      = ncol(df)
    @assert ncols ≥ 5 "Expect at least :iter plus one 4-column FW block"
    n_variants = Int((ncols - 1) ÷ 4)

    # ---------------- build human-readable labels (unchanged) ----------------------
    FW_TOKENS      = ["AP","BPFW","AFW","FW","BC","C","F"]
    token_re       = Regex("($(join(sort(FW_TOKENS, by=length, rev=true), '|')))")

    variant_labels = String[]
    for v in 1:n_variants
        prim_header = String(cols[1 + (v-1)*4 + 1])  # after `iter`, offset +1
        block_name  = split(prim_header, "_")[1]
        push!(variant_labels, join([m.match for m in eachmatch(token_re, block_name)], "-"))
    end

    # ---------------- keep only requested variants ---------------------------------
    wanted_idxs = isempty(wanted_fw_variants) ?
                  collect(1:n_variants) :
                  findall(lbl -> lbl in wanted_fw_variants, variant_labels)
    @assert !isempty(wanted_idxs) "No matching FW variants: $wanted_fw_variants"

    labels       = variant_labels[wanted_idxs]
    trajectories = [Vector{Any}(undef, nrow(df)) for _ in wanted_idxs]

    # ---------------- fill trajectories --------------------------------------------
    for row_idx in 1:nrow(df)
        iter = Int64(df[row_idx, cols[1]])

        for local_idx in 1:length(wanted_idxs)
            orig_idx = wanted_idxs[local_idx]
            base     = 1 + (orig_idx - 1)*4

            prim = Float64(df[row_idx, cols[base + 1]])
            dual = Float64(df[row_idx, cols[base + 2]])
            dgap = Float64(df[row_idx, cols[base + 3]])
            time = Float64(df[row_idx, cols[base + 4]])

            trajectories[local_idx][row_idx] = (iter, prim, dual, dgap, time)
        end
    end

    return trajectories, labels 
end





"""
    average_over_logs(log_trajectories; fw_variants_labels::Vector{String})
Input:
    - `log_trajectories`: contains all different instances, each with a traj element. Each traj is Vector{Vector{Tuple()}}
    - `fw_variants_labels`: labels of FW variants analyzed in each trajectory
Output: 
    - `averaged_trajectories::Vector{Vector{Any}}`: Vector containing, ∀ FW variant, a Vector of tuples avged across all logfiles,
       sorted in ascending `iter`
"""
function avg_over_logs(log_trajectories::Vector{Vector{Vector{Any}}}, fw_variants_labels::Vector{String})
    m = length(fw_variants_labels)
    L = length(log_trajectories)

    # --- 1) Build a DataFrame for each logfile j ---
    dfs = Vector{DataFrame}(undef, L)
    for j in 1:L
        trajs = log_trajectories[j]  # Vector of length m

        # Prepare column‐vectors
        logfile = Int[]
        variant = String[]
        iter    = Int[]
        pgap    = Float64[]
        dual    = Float64[]
        dgap    = Float64[]
        time    = Float64[]

        for i in 1:m
            lbl = fw_variants_labels[i]
            for (it, pg, du, dg, tm) in trajs[i]
                push!(logfile, j)
                push!(variant, lbl)
                push!(iter,    it)
                push!(pgap,    pg)
                push!(dual,    du)
                push!(dgap,    dg)
                push!(time,    tm)
            end
        end

        dfs[j] = DataFrame(
            logfile = logfile,
            variant = variant,
            iter    = iter,
            pgap    = pgap,
            dual    = dual,
            dgap    = dgap,
            time    = time,
        )
    end

    # --- 2) Concatenate and group‐average ---
    bigdf = vcat(dfs...)
    g     = groupby(bigdf, [:variant, :iter])
    avgdf = combine(g,
        :pgap => mean => :pgap_mean,
        :dual => mean => :dual_mean,
        :dgap => mean => :dgap_mean,
        :time => mean => :time_mean
    )
    # avgdf columns: :variant, :iter, :pgap_mean, :dual_mean, :dgap_mean, :time_mean

    # --- 3) Reconstruct into Vector{Vector{Tuple}} order by fw_variants_labels ---
    averaged_trajectories = Vector{Vector{Tuple{Int,Float64,Float64,Float64,Float64}}}(undef, m)
    for (i, lbl) in enumerate(fw_variants_labels)
        sub = @view avgdf[avgdf.variant .== lbl, :]
        sort!(sub, :iter)  # ensure increasing iteration
        averaged_trajectories[i] = [(row.iter,
                                     row.pgap_mean,
                                     row.dual_mean,
                                     row.dgap_mean,
                                     row.time_mean) for row in eachrow(sub)]
    end

    return averaged_trajectories
end