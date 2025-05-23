# plot_from_logs.jl
# This example shows how to build plots from an existing logfile. 
# This is useful if the user decides at a certain point to change anything in the plotting functions. but does not want
#   At that point, instead of re-running all the experiments from scratch, one can just plot using saved logs
#   One can decide to plot only some FW variants from the logs
using ProductPolytopesAFW
using FrankWolfe
using Plots



# ---------------------------------------------------------------------------------
# YAML PARAMETERS
# ---------------------------------------------------------------------------------
config = Config("examples/config.yml")

# ---------------------------------------------------------------------------------
# SCRIPT PARAMETERS
# ---------------------------------------------------------------------------------
basename_log = "i_k2_n10000_i1000_s240389_cvxho_anc_t20250511_123700"
logname = basename_log*".csv"

k, n = get_k_n_from_logstring(basename_log)
config = modify_config(config, k=k, n=n)

# ---------------------------------------------------------------------------------
# RETRIEVE AND PROCESS LOGS
# ---------------------------------------------------------------------------------
# retrieve logs
wanted = ["C-BC-FW", "F-FW", "F-AFW"]

# ---------------------------------------------------------------------------------
# PLOTS NON-INTERSECTING INSTANCES
# trajis_ni, labels, opt = load_fw_trajectories_ni("examples/results_linesearch_afw/iter_logs/$logname", wanted_fw_variants=wanted)
# # find `cutoff_time` of the FW variant that ends first, then cutoff all iters of all other variants happening after `cutoff_time`
# cutoff_trajectories_ni, _ = cutoff_log_shortest_time(trajis_ni)
# # compute primal gap from primal, for all nonintersecting instances
# cutoff_trajectories_ni_pgap = [compute_primal_gap(t, opt) for t in cutoff_trajectories_ni]
# # plot only primal and FW gap over time
# fig_ni = plot_time_only(config, cutoff_trajectories_ni_pgap, labels, yscalelog=true, xscalelog=true)
# # decide figure name
# fig_ni_filename = "examples/results_linesearch_afw/plots/plot_$(basename_log)_fromlogs.pdf"
# # Plot trajectories and save as PDF
# Plots.savefig(fig_ni, fig_ni_filename)

# ---------------------------------------------------------------------------------
# PLOTS INTERSECTING INSTANCES
trajis_i, labels = load_fw_trajectories_i("examples/results_linesearch_afw/iter_logs/$logname", wanted_fw_variants=wanted)
cutoff_trajectories_i, _ = cutoff_log_shortest_time(trajis_i)
fig_i  = plot_time_only(config, cutoff_trajectories_i, labels, yscalelog=true, xscalelog=true)
fig_i_filename = "examples/results_linesearch_afw/plots/plot_$(basename_log)_fromlogs.pdf"
Plots.savefig(fig_i, fig_i_filename)