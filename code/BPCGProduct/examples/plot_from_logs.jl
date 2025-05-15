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
basename = "ni_k2_n20000_i1000_s240389_cvxho_anc_t20250511_212230"
logname = basename*".csv"

k, n = get_k_n_from_logstring(basename)
config = modify_config(config, k=k, n=n)

# ---------------------------------------------------------------------------------
# RETRIEVE AND PROCESS LOGS
# ---------------------------------------------------------------------------------
# retrieve logs
wanted = ["C-BC-FW", "F-BC-AFW", "F-FW", "F-AFW"]
trajis_ni, labels, opt = load_fw_trajectories("examples/results_linesearch_afw/iter_logs/$logname", wanted_fw_variants=wanted)
println(labels)
# find `cutoff_time` of the FW variant that ends first, then cutoff all iters of all other variants happening after `cutoff_time`
cutoff_trajectories_ni, cutoff_time_ni = cutoff_log_shortest_time(trajis_ni)
# compute primal gap from primal, for all nonintersecting instances
cutoff_trajectories_ni_pgap = [compute_primal_gap(t, opt) for t in cutoff_trajectories_ni]

# ---------------------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------------------
# plot only primal and FW gap over time
fig_ni = plot_time_only(config, cutoff_trajectories_ni_pgap, labels, yscalelog=true, xscalelog=true)
# decide figure name
fig_ni_filename = "examples/results_linesearch_afw/plots/plot_$(basename)_fromlogs.pdf"
# Plot trajectories and save as PDF
Plots.savefig(fig_ni, fig_ni_filename)