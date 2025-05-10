# plot_from_logs.j
# This example shows how to build plots from existing plots. 
# This is useful if the user decides at a certain point to change anything in the plotting functions. but does not want
#   At that point, instead of re-running all the experiments from scratch, one can just plot using saved logs
using BPCGProduct
using FrankWolfe
using Plots

 
basename = "ni_k3_n101_s15672_cvxho_anc_t20250423170234"
logname = basename*".csv"

# retrieve logs
trajis_ni, variant_labels, opt = load_fw_trajectories("examples/results_linesearch_afw/logs/$logname")
# find `cutoff_time` of the FW variant that ends first, then cutoff all iters of all other variants happening after `cutoff_time`
cutoff_trajectories_ni, cutoff_time_ni = cutoff_log_shortest_time(trajis_ni)
# plot only primal and FW gap over time
fig_ni = plot_time_only(cutoff_trajectories_ni, variant_labels, yscalelog=true, xscalelog=true)
# decide figure name
fig_ni_filename = "examples/results_linesearch_afw/plots/plot_$basename.pdf"
# Plot trajectories and save as PDF
Plots.savefig(fig_ni, fig_ni_filename)