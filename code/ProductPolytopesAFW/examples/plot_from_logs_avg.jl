# plot_from_logs_avg.jl
# This example shows how to build plots from existing logfiles, averaging their values
# This is useful if the user decides at a certain point to change anything in the plotting functions. but does not want
#   At that point, instead of re-running all the experiments from scratch, one can just plot using saved logs
#   One can decide to plot only some FW variants from the logs
using ProductPolytopesAFW
using FrankWolfe
using Plots


"""
    files_with_prefix(dir::AbstractString, prefix::AbstractString; fullpaths::Bool=false)

Return a Vector of filenames in `dir` whose names start with `prefix` 
"""
function retrieve_logfiles_by_prefix(dir::AbstractString, prefix::AbstractString)
    
    # `readdir(dir)`: read all entries (files and subdirectories) from `dir`, then
    # filter those a) starting w/given prefix, b) ending with ".csv" and c) being files (no directories)
    # finally, `joinpath` gives the full path for each of the objects
    return [
        joinpath(dir, f) for f in readdir(dir) 
            if startswith(f, prefix) && endswith(f, ".csv") && isfile(joinpath(dir, f))
        ]
end

using Statistics

"""
    average_fw_trajectories_ni(logdir::AbstractString,
                              prefix::AbstractString,
                              wanted_fw_variants::Vector{String})

Find every “ni_*…*.csv” in `logdir` whose filename starts with `prefix`, then:
 1. Load them with `load_fw_trajectories_ni`.
 2. Compute the per-log cutoff time (earliest-finishing variant).
 3. Find the global minimum of those cutoff times.
 4. Re-truncate *every* trajectory in every log at that global cutoff.
 5. Compute primal-gap on each truncated trajectory.
 6. For each FW variant, average the (iter, pgap, dual, dgap, time) tuples across logs.

Returns `(avg_trajs, labels)` where
- `avg_trajs` is a `Vector` of length `length(wanted_fw_variants)`, each an
  `Vector{NTuple{5,Float64}}` of the averaged tuples;
- `labels` is the human-readable names of the variants.
"""
function average_fw_trajectories_ni(logdir::AbstractString,
                                    prefix::AbstractString,
                                    wanted_fw_variants::Vector{String})

    # 1) gather the filenames
    files = files_with_prefix(logdir, prefix)
    @assert !isempty(files) "No log files found matching “$prefix” in $logdir"

    # 2) for each file, load and compute that log’s cutoff time
    per_log_cutoffs = Float64[]
    for file in files
        trajs, labels, opt = load_fw_trajectories_ni(file;
                                  wanted_fw_variants=wanted_fw_variants)
        _, cutoff = cutoff_log_shortest_time(trajs)
        push!(per_log_cutoffs, cutoff)
    end

    # 3) global cutoff time
    global_cutoff = minimum(per_log_cutoffs)

    # 4) re‐truncate each log at global_cutoff and compute primal gap
    #    collect per‐log “pgap” trajectories into a 3D array:
    #      dims = (n_variants, n_iters, n_logs)
    n_variants = length(wanted_fw_variants)
    # load one to figure out how many iterations after truncation:
    example_trajs, _, example_opt = load_fw_trajectories_ni(files[1];
                                         wanted_fw_variants=wanted_fw_variants)
    truncated_example, _ = cutoff_log_shortest_time(example_trajs,
                                                    global_cutoff)
    n_iters = length(truncated_example[1])

    # space to accumulate
    pgap_arr   = zeros(n_variants, n_iters, length(files))
    dual_arr   = zeros(n_variants, n_iters, length(files))
    dgap_arr   = zeros(n_variants, n_iters, length(files))
    time_arr   = zeros(n_variants, n_iters, length(files))
    iters      = Vector{Int}(undef, n_iters)  # iteration indices

    # fill them
    for (i_log, file) in ipairs(files)
        trajs, _, opt = load_fw_trajectories_ni(file;
                                  wanted_fw_variants=wanted_fw_variants)
        trunc, _ = cutoff_log_shortest_time(trajs, global_cutoff)
        for v in 1:n_variants, t in 1:n_iters
            iter, primal, dual, dgap, time = trunc[v][t]
            iters[t]      = iter
            pgap_arr[v,t,i_log] = primal - opt
            dual_arr[v,t,i_log] = dual
            dgap_arr[v,t,i_log] = dgap
            time_arr[v,t,i_log] = time
        end
    end

    # 5) average across the log‐dimension
    avg_trajs = Vector{Vector{NTuple{5,Float64}}}(undef, n_variants)
    for v in 1:n_variants
        avg_trajs[v] = NTuple{5,Float64}[]
        for t in 1:n_iters
            push!(avg_trajs[v], (
                iters[t],
                mean(pgap_arr[v,t,:]),
                mean(dual_arr[v,t,:]),
                mean(dgap_arr[v,t,:]),
                mean(time_arr[v,t,:]),
            ))
        end
    end

    return avg_trajs, wanted_fw_variants
end




# ---------------------------------------------------------------------------------
# YAML PARAMETERS
# ---------------------------------------------------------------------------------
config = Config("examples/config.yml")

# ---------------------------------------------------------------------------------
# SCRIPT PARAMETERS
# ---------------------------------------------------------------------------------
lognames = retrieve_logfiles_by_prefix("examples/results_linesearch_afw/iter_logs", "ni_k2_n10000_i1000_s")
basename = "ni_k2_n20000_i1000_cvxho_anc_avg$(length(lognames))s"
k, n = get_k_n_from_logstring(basename)
config = modify_config(config, k=k, n=n)

# ---------------------------------------------------------------------------------
# RETRIEVE AND PROCESS LOGS
# ---------------------------------------------------------------------------------
# retrieve logs
wanted = ["C-BC-FW", "F-FW", "F-AFW"]