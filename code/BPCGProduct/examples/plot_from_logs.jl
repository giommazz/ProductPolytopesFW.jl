# plot_from_logs.j
# This example shows how to build plots from existing plots. 
# This is useful if the user decides at a certain point to change anything in the plotting functions. but does not want
#   At that point, instead of re-running all the experiments from scratch, one can just plot using saved logs
using BPCGProduct
using FrankWolfe
using Plots

# ---------------------------------------------------------------------------------
# YAML PARAMETERS
# ---------------------------------------------------------------------------------
config = Config("examples/config.yml")

# ---------------------------------------------------------------------------------
# SCRIPT PARAMETERS
# ---------------------------------------------------------------------------------
basename = "i_k2_n20000_i1000_s240389_cvxho_anc_t20250511_212230"# "ni_k5_n10000_i1000_s240389_cvxho_anc_t20250512_140502"# #"ni_k2_n20000_i1000_s240389_cvxho_anc_t20250511_212230"#, "ni_k10_n10000_i1000_s240389_cvxho_anc_t20250513_034433"
logname = basename*".csv"

k, n = get_k_n_from_logstring(basename)
config = modify_config(config, k=k, n=n)

# ---------------------------------------------------------------------------------
# RETRIEVE AND PROCESS LOGS
# ---------------------------------------------------------------------------------
# retrieve logs
wanted = ["C-BC-FW", "F-FW", "F-AFW"]
"""
iter,CBCFW_prim,CBCFW_dual,CBCFW_dgap,CBCFW_time,FBCFW_prim,FBCFW_dual,FBCFW_dgap,FBCFW_time,FBCAFW_prim,FBCAFW_dual,FBCAFW_dgap,FBCAFW_time,FFW_prim,FFW_dual,FFW_dgap,FFW_time,FAFW_prim,FAFW_dual,FAFW_dgap,FAFW_time
"""

# *********************************************************************************
# *********************************************************************************
# *********************************************************************************
# *********************************************************************************
# *********************************************************************************
# *********************************************************************************
using CSV, DataFrames
function loady_ni(path::String; wanted_fw_variants::Vector{String}=String[])
    
    # `splitext` returns (root, extension)
    ext = lowercase(splitext(path)[2])
    @assert ext == ".csv" "Expected a .csv file"
    
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
        
        # for each variant, grab its 4 fields
        # for FWvar_index in 1:length(retrieved_fw_labels)
        #     # skip the first two columns (iter,opt), then block of 4 per FW variant
        #     base = 2 + (FWvar_index - 1)*4
        #     # ...then retrieve tuple indices for each variant and each row
        #     prim_col = cols[base+1]
        #     dual_col = cols[base+2]
        #     dgap_col = cols[base+3]
        #     time_col = cols[base+4]
        #     # now you have row indices and column indices, convert values from String to Number
        #     prim = Float64(df[row_idx, prim_col])
        #     dual = Float64(df[row_idx, dual_col])
        #     dgap = Float64(df[row_idx, dgap_col])
        #     time = Float64(df[row_idx, time_col])

        #     trajectories[FWvar_index][row_idx] = (iter, prim, dual, dgap, time)
        # end

        for local_idx in 1:length(wanted_idxs)
            orig_idx = wanted_idxs[local_idx]
            base     = 2 + (orig_idx - 1) * 4

            prim = Float64(df[row_idx, cols[base + 1]])
            dual = Float64(df[row_idx, cols[base + 2]])
            dgap = Float64(df[row_idx, cols[base + 3]])
            time = Float64(df[row_idx, cols[base + 4]])

            trajectories[local_idx][row_idx] = (iter, prim, dual, dgap, time)
        end
    end

    return trajectories, retrieved_fw_labels, opt
end





function loady_i(path::String; wanted_fw_variants::Vector{String}=String[])
    # ---------------- basic checks -------------------------------------------------
    @assert lowercase(splitext(path)[2]) == ".csv" "Expected a .csv file"

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



function loady_fw_trajectories(path::String; wanted_fw_variants::Vector{String}=String[])
    
    # `splitext` returns (root, extension)
	@assert lowercase(splitext(path)[2]) == ".csv" "Expected a .csv file"

    # read .csv file
    df = CSV.read(path, DataFrame)
    # extract column names
    cols = names(df)
    # determine n. of FW variants that were executed (columns are :iter,:opt, then 4×n_variants columns)
    ncols = ncol(df)
	
	@assert (length(ncols) % 4 == 1 || length(ncols) % 4 == 2) "there is something wrong with the columns in $path"
	if "opt" in names(df)	# non-intersecting instance, "opt" column present
		opt_col_counter = 1
	else 					# intersecting instance, no "opt" column
		opt_col_counter = 0
	end
	
    expected_ncols = 5 + opt_col_counter
	@assert ncols ≥ expected_ncols "Expected ≥ $expected_ncols columns (iter, opt, plus at least one variant's 4 columns)"
    n_variants = Int((ncols - 1 - opt_col_counter) ÷ 4)

    
    # define allowed tokens
    FW_TOKENS = ["AP", "BPFW", "AFW", "FW", "BC", "C", "F"]
    # sort tokens by descending length so that multi-letter tokens (e.g. "BPFW") match before shorter ones ("FW")
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
        prim_header = String(cols[1 + opt_col_counter + (v-1)*4 + 1])         # e.g. "CBCFW_prim"
        block_name  = split(prim_header, "_")[1]    # e.g. "CBCFW"

        # d) find all occurrences of tokens in that block name, yielding ["C","BC","FW"]
        parts = [m.match for m in eachmatch(token_re, block_name)]
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
    opt = opt_col_counter > 0 ? Float64(df[1, cols[2]]) : 0.0

    # pre-allocate the output: one vector per variant
    trajectories = Vector{Vector{Any}}(undef, length(retrieved_fw_labels))
    for i in 1:length(retrieved_fw_labels)
        trajectories[i] = Vector{Any}(undef, nrow(df))
    end

    # row indices
    for row_idx in 1:nrow(df)
        
        # read iteration number
        iter = Int64(df[row_idx, cols[1]])
        
        # for each variant, grab its 4 fields
        # for FWvar_index in 1:length(retrieved_fw_labels)
        #     # skip the first two columns (iter,opt), then block of 4 per FW variant
        #     base = 2 + (FWvar_index - 1)*4
        #     # ...then retrieve tuple indices for each variant and each row
        #     prim_col = cols[base+1]
        #     dual_col = cols[base+2]
        #     dgap_col = cols[base+3]
        #     time_col = cols[base+4]
        #     # now you have row indices and column indices, convert values from String to Number
        #     prim = Float64(df[row_idx, prim_col])
        #     dual = Float64(df[row_idx, dual_col])
        #     dgap = Float64(df[row_idx, dgap_col])
        #     time = Float64(df[row_idx, time_col])

        #     trajectories[FWvar_index][row_idx] = (iter, prim, dual, dgap, time)
        # end

        for local_idx in 1:length(wanted_idxs)
            orig_idx = wanted_idxs[local_idx]
            base     = 1 + opt_col_counter + (orig_idx - 1) * 4

            prim = Float64(df[row_idx, cols[base + 1]])
            dual = Float64(df[row_idx, cols[base + 2]])
            dgap = Float64(df[row_idx, cols[base + 3]])
            time = Float64(df[row_idx, cols[base + 4]])

            trajectories[local_idx][row_idx] = (iter, prim, dual, dgap, time)
        end
    end

    return trajectories, retrieved_fw_labels, opt
end



# PLOTS NI
# ---------------------------------------------------------------------------------
# trajis_ni, labels, opt = loady_ni("examples/results_linesearch_afw/iter_logs/$logname", wanted_fw_variants=wanted)
# # find `cutoff_time` of the FW variant that ends first, then cutoff all iters of all other variants happening after `cutoff_time`
# cutoff_trajectories_ni, _ = cutoff_log_shortest_time(trajis_ni)
# # compute primal gap from primal, for all nonintersecting instances
# cutoff_trajectories_ni_pgap = [compute_primal_gap(t, opt) for t in cutoff_trajectories_ni]
# # plot only primal and FW gap over time
# fig_ni = plot_time_only(config, cutoff_trajectories_ni_pgap, labels, yscalelog=true, xscalelog=true)
# # decide figure name
# fig_ni_filename = "examples/results_linesearch_afw/plots/plot_$(basename)_fromlogs.pdf"
# # Plot trajectories and save as PDF
# Plots.savefig(fig_ni, fig_ni_filename)


# PLOTS NI
# ---------------------------------------------------------------------------------
trajis_i, labels = loady_i("examples/results_linesearch_afw/iter_logs/$logname", wanted_fw_variants=wanted)
println(labels)
readline()
cutoff_trajectories_i, _ = cutoff_log_shortest_time(trajis_i)
fig_i  = plot_time_only(config, cutoff_trajectories_i, labels, yscalelog=true, xscalelog=true)
fig_i_filename = "examples/results_linesearch_afw/plots/plot_$(basename)_fromlogs.pdf"
Plots.savefig(fig_i, fig_i_filename)