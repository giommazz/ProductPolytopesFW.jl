# plot_from_logs_avg.jl
# This example shows how to build plots from existing logfiles, averaging their values
# This is useful if the user decides at a certain point to change anything in the plotting functions. but does not want
#   At that point, instead of re-running all the experiments from scratch, one can just plot using saved logs
#   One can decide to plot only some FW variants from the logs
using ProductPolytopesAFW
using FrankWolfe
using Plots


# --------------------------------------------------------------
# HELPER FUNCTIONS
# --------------------------------------------------------------
"""
    retrieve_logfiles_by_prefix(dir::String, prefix::String)
Input: 
    - logdir: directory containing the logfiles
    - prefix: string used to search for specific logfiles

Output: 
    - `logfiles`: Vector of filenames in `logdir` whose names start with `prefix` 
"""
function retrieve_logfiles_by_prefix(logdir::String, prefix::String)
    
    # `readdir(dir)`: read all entries (files and subdirectories) from `dir`, then
    # filter those a) starting w/given prefix, b) ending with ".csv" and c) being files (no directories)
    # finally, `joinpath` gives the full path for each of the objects
    logfiles =  [
        joinpath(logdir, f) for f in readdir(logdir)
            if startswith(f, prefix) && endswith(f, ".csv") && isfile(joinpath(logdir, f))
        ]
    @assert !isempty(logfiles) "No log files found matching $prefix in $logdir"
    return logfiles
end


"""
    average_fw_trajectories(logdir::AbstractString, prefix::AbstractString, wanted_fw_variants::Vector{String})
"""
function average_fw_trajectories(
    logfiles::Vector{String},
    prefix::String,
    wanted_fw_variants::Vector{String},
    ni_flag::Bool)

    # compute global cutoff time, i.e., min cutoff time over all trajectories in all logfiles
    global_cutoff_time = cutoff_time(logfiles, wanted_fw_variants, ni_flag)    

    # load, then get the the trajectories list and the optimal value list (if non-intersecting)
    if ni_flag
        loaded_results = [load_fw_trajectories_ni(logfile; wanted_fw_variants=wanted_fw_variants) for logfile in logfiles]
        opts_all_logs = getindex.(loaded_results, 3)
    else
        loaded_results = [load_fw_trajectories_i(logfile; wanted_fw_variants=wanted_fw_variants) for logfile in logfiles]
    end

    traj_all_logs = getindex.(loaded_results, 1)

    # compute cutoff trajectories based on `global_cutoff_time`
    cutoff_trajectories_all_logs = [cutoff_log_shortest_time(traj, global_cutoff_time) for traj in traj_all_logs]

    if ni_flag
        cutoff_trajectories_pgap_all_logs = [
            # contains all updated pgaps for all of the `wanted_fw_variants` in this logfile's instance
            [compute_primal_gap(cutoff_trajectories_all_logs[logfile_idx][fw_var], opts_all_logs[logfile_idx]) 
                # for each FW variant
                for fw_var in 1:length(wanted_fw_variants)]
            # for each instance in the logfiles
            for logfile_idx in 1:length(traj_all_logs)
            ]
    else
        cutoff_trajectories_pgap_all_logs = cutoff_trajectories_all_logs
    end

    avg_cutoff_trajectories_pgap = avg_over_logs(cutoff_trajectories_pgap_all_logs, wanted_fw_variants)

    # plot only primal and FW gap over time
    figg = plot_time_only(config, avg_cutoff_trajectories_pgap, wanted_fw_variants, yscalelog=true, xscalelog=true)
    # decide figure name
    figg_filename = "examples/results_linesearch_afw/plots/plot_$(basename)avg.pdf"
    # Plot trajectories and save as PDF
    Plots.savefig(figg, figg_filename)

    return avg_cutoff_trajectories_pgap
end





# ---------------------------------------------------------------------------------
# YAML PARAMETERS
# ---------------------------------------------------------------------------------
config = Config("examples/config.yml")


# ---------------------------------------------------------------------------------
# SCRIPT PARAMETERS
# ---------------------------------------------------------------------------------
prefix = "ni_k2_n20000_i1000"
ni_flag = prefix[1:3] == "ni_" ? true : false

logdir = "examples/results_linesearch_afw/iter_logs"
wanted_fw_variants = ["C-BC-FW", "F-FW", "F-AFW"]

# gather the filenames
logfiles = retrieve_logfiles_by_prefix(logdir, prefix)
# decide basename and modify config (needed for plot legends)
basename_avg_plot = prefix*"_cvxho_anc_avg$(length(logfiles))"
k, n = get_k_n_from_logstring(basename_avg_plot)
config = modify_config(config, k=k, n=n)

# retrieve and process trajectories, then compute averages
avg_cutoff_trajectories_pgap = average_fw_trajectories(logfiles, prefix, wanted_fw_variants, ni_flag)

# plot only primal and FW gap over time
figg = plot_time_only(config, avg_cutoff_trajectories_pgap, wanted_fw_variants, yscalelog=true, xscalelog=true)
# decide figure name
figg_filename = "examples/results_linesearch_afw/plots/plot_$(basename_avg_plot).pdf"
# Plot trajectories and save as PDF
Plots.savefig(figg, figg_filename)