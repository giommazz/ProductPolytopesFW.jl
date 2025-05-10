# `compute_intersection_custom_full_warmup_slurm.jl`
# Run either within the Julia REPL as include("/examples/compute_intersection_custom_full.jl")
# or from linux terminal with: `julia --project=. examples/compute_intersection_custom_full.jl > test.log 2>&1`
#
# This includes a little script that warms up the REPL: it is executed before running the main script, so that when the main script is run, 
#   the package `BPCGProduct` has already been precompiled and the plots obtained by `main` don't show initial arbitrary overhead
using BPCGProduct
using FrankWolfe

# ---------------------------------------------------------------------------------
# MAIN FUNCTIONS
# Solve small instance (using `config_warmup`) to "warm-up" the REPL: this compiles `BPCGProduct`, so that no compilation needed upon running `main`
function repl_warmup(config::Config, vertices, shifted_vertices, labels, basename)

    # Retrieve nonintersecting and intersecting LMOs from previously generated instances
    lmo_list = create_lmos(config, [vertices, shifted_vertices])

    # Will contain data about diafferent FW runs, for non-intersecting and intersecting polytopes
    trajectories_ni, trajectories_i = [], []

    for (i, lmos) in enumerate(lmo_list)
        # nonintersecting flag
        ni_flag = i == 1
        println()
        if ni_flag println("\nNon intersecting") else println("\nIntersecting") end
        prod_lmo = create_product_lmo(lmos)
        _, _, _, _, _ = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        _, _, _, _, _ = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        _, _, _, _, _ = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), AwayStep(), prod_lmo)  
        _, _, _, _, _ = run_FullFW(config, FrankWolfe.frank_wolfe, prod_lmo)    
        _, _, _, _, _ = run_FullFW(config, FrankWolfe.away_frank_wolfe, prod_lmo)    
    end

    return trajectories_ni, trajectories_i
end

function main(config::Config, vertices, shifted_vertices, opt, labels, basename)

    # Retrieve nonintersecting and intersecting LMOs from previously generated instances
    lmo_list = create_lmos(config, [vertices, shifted_vertices])

    # Will contain data about diafferent FW runs, for non-intersecting and intersecting polytopes
    trajectories_ni, trajectories_i = [], []

    for (i, lmos) in enumerate(lmo_list)
        # nonintersecting flag
        ni_flag = i == 1
        println()
        if ni_flag println("\nNon intersecting") else println("\nIntersecting") end

        prod_lmo = create_product_lmo(lmos)
        
        # Run Frank-Wolfe algorithms and alternating projections, then record trajectory data
        
        # ***************************************
        # Cyclic block-coordinate methods
        println("\n\n\n ----------> Cyclic Block-coordinate vanilla FW")
        _, _, _, _, td_cyc_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_cyc_bc_fw, trajectories_ni, trajectories_i, opt)
        
        # println("\n\n\n ----------> Cyclic Block-coordinate AFW")
        # _, _, _, _, td_cyc_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), AwayStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_cyc_bc_afw, trajectories_ni, trajectories_i, opt)

        # println("\n\n\n ----------> Stochastic Block-coordinate vanilla FW")
        # _, _, _, _, td_stoc_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.StochasticUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_stoc_bc_fw, trajectories_ni, trajectories_i, opt)
        
        # println("\n\n\n ----------> Stochastic Block-coordinate AFW")
        # _, _, _, _, td_stoc_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.StochasticUpdate(), AwayStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_stoc_bc_afw, trajectories_ni, trajectories_i, opt)
        
        # println("\n\n\n ----------> Cyclic Block-coordinate BPFW")
        # _, _, _, _, td_cyc_bc_bpfw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_cyc_bc_bpfw, trajectories_ni, trajectories_i, opt)

        # ***************************************
        # Full block-coordinate methods
        println("\n\n\n ----------> Full Block-coordinate vanilla FW")
        _, _, _, _, td_full_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_full_bc_fw, trajectories_ni, trajectories_i, opt)

        println("\n\n\n ----------> Full Block-coordinate Away FW (ours)")
        _, _, _, _, td_full_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), AwayStep(), prod_lmo)  
        push_to_trajectories!(ni_flag, td_full_bc_afw, trajectories_ni, trajectories_i, opt)

        # println("\n\n\n ----------> Full Block-coordinate Blended Pairwise FW (ours)")
        # _, _, _, _, td_full_bc_bpfw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), FrankWolfe.BPCGStep(), prod_lmo)  
        # push_to_trajectories!(ni_flag, td_full_bc_bpfw, trajectories_ni, trajectories_i, opt)

        # ***************************************
        # Full methods
        println("\n\n\n ----------> Full FW")
        _, _, _, _, td_full_fw = run_FullFW(config, FrankWolfe.frank_wolfe, prod_lmo)    
        push_to_trajectories!(ni_flag, td_full_fw, trajectories_ni, trajectories_i, opt)

        println("\n\n\n ----------> Full AFW")
        _, _, _, _, td_full_afw = run_FullFW(config, FrankWolfe.away_frank_wolfe, prod_lmo)    
        push_to_trajectories!(ni_flag, td_full_afw, trajectories_ni, trajectories_i, opt)

        # println("\n\n\n ----------> Full BPFW")
        # _, _, _, _, td_full_bpfw = run_FullFW(config, FrankWolfe.blended_pairwise_conditional_gradient, prod_lmo)    
        # push_to_trajectories!(ni_flag, td_full_bpfw, trajectories_ni, trajectories_i, opt)
        
        # println("\n\n\n ----------> AP")
        # _, _, _, _, td_ap = run_AlternatingProjections(config, prod_lmo, true)    
        # # `FrankWolfe.alternating_projections` computes ||x-y|| rather than 1/2 ||x-y||
        # push_to_trajectories!(ni_flag, td_ap, trajectories_ni, trajectories_i, 2*opt)
    end

    return trajectories_ni, trajectories_i
end





t_start_tot = time()





# ---------------------------------------------------------------------------------
# YAML PARAMETERS
# ---------------------------------------------------------------------------------
config = Config("examples/config.yml")
print_config(config)
config_warmup = modify_config(config, k=2, n=15)
print_config(config_warmup)





# ---------------------------------------------------------------------------------
# WARM-UP SCRIPT
# ---------------------------------------------------------------------------------
# Generate instances 
println()
println()
println("********************************************************")
println("WARMUP: Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, _, fw_gap = generate_polytopes(config_warmup)
# Optimal solution
basename = generate_filename(config_warmup)
# Labels for the plots
labels = ["F-BC-AFW"]
# execute main
println("********************************************************")
println("WARMUP: Running FW on the instances")
println("********************************************************")
_, _ = repl_warmup(config_warmup, vertices, shifted_vertices, labels, basename)





# ---------------------------------------------------------------------------------
# MAIN SCRIPT
# ---------------------------------------------------------------------------------
println()
println()
println()
# ---------------------------------------------------------------------------------
# PARAMETERS
results_dir = "examples/results_linesearch_afw" # "results_shortstep", 
times_dir = results_dir*"/times"
logs_dir = results_dir*"/iter_logs"
plots_dir = results_dir*"/plots"
# Labels for the plots
labels = ["C-BC-FW", "F-BC-FW", "F-BC-AFW", "F-FW", "F-AFW"] # ["C-BC-FW", "C-BC-AFW", "C-BC-BPFW", "F-BC-FW", "F-BC-AFW", "F-BC-BPFW", "F-FW", "F-AFW", "F-BPFW", "AP"]

# ---------------------------------------------------------------------------------
# GENERATE INSTANCES
println("********************************************************")
println("MAIN: Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, opt, fw_gap = generate_polytopes(config)
basename = generate_filename(config)

# ---------------------------------------------------------------------------------
# RUN MAIN TO SOLVE INSTANCES
println("\n\n********************************************************")
println("MAIN: Running FW on the instances")
println("********************************************************")
t_start_main = time()
trajectories_ni, trajectories_i = main(config, vertices, shifted_vertices, opt, labels, "ni_"*basename)
t_end_main = time()
println("\n\n\t\tElapsed time main: ", t_end_main - t_start_main, " seconds\n\n")

# ---------------------------------------------------------------------------------
# PROCESS AND SAVE LOGS
# Save log data in `.csv` format
# pad data so that all FW runs have the same number of iterations/lines
padded_trajectories_ni, min_length_ni, max_length_ni = pad_log_data(trajectories_ni)
padded_trajectories_i, min_length_i, max_length_i = pad_log_data(trajectories_i)
# log padded data
save_logdata_to_csv(padded_trajectories_ni, opt, max_length_ni, labels, logs_dir, "ni_"*basename)
save_logdata_to_csv(padded_trajectories_i, max_length_i, labels, logs_dir, "i_"*basename)

# ---------------------------------------------------------------------------------
# CREATE AND SAVE PLOTS
# cutoff trajectories at the same iteration (shortest run determines the iteration), for plotting
cutoff_trajectories_ni, _ = cutoff_log_shortest_time(padded_trajectories_ni)
cutoff_trajectories_i, _ = cutoff_log_shortest_time(padded_trajectories_i)

# compute primal gap from primal, for all nonintersecting instances
cutoff_trajectories_ni_pgap = [compute_primal_gap(t, opt) for t in cutoff_trajectories_ni]


# Generate plots (do not pass `filename` argument, so .png is not automatically saved)
# fig_ni = plot_trajectories(trajectories_ni, labels, yscalelog=true, xscalelog=true) # plotting function from FrankWolfe.jl package
# fig_i = plot_trajectories(trajectories_i, labels, yscalelog=true, xscalelog=true) # plotting function from FrankWolfe.jl package
fig_ni = plot_time_only(config, cutoff_trajectories_ni_pgap, labels, yscalelog=true, xscalelog=true) # plotting function that only prints time and is customized
fig_i  = plot_time_only(config, cutoff_trajectories_i, labels, yscalelog=true, xscalelog=true) # plotting function that only prints time and is customized
# Decide filename
fig_ni_filename = plots_dir*"/plot_ni_$basename"
fig_i_filename = plots_dir*"/plot_i_$basename"
# Plot trajectories
# Plots.plot!(fig_ni, size=(1200, 800))  # Larger figure size
# Plots.plot!(fig_i, size=(1200, 800))  # Larger figure size
# Manually save only the PDF versions
Plots.savefig(fig_ni, fig_ni_filename*".pdf")  
Plots.savefig(fig_i, fig_i_filename*".pdf")

# ---------------------------------------------------------------------------------
# SAVE TIMES
# Save time data in `.csv` format
log_times(trajectories_i, labels, times_dir*"/times_i_"*basename)
log_times(trajectories_ni, labels, times_dir*"/times_ni_"*basename)





t_end_tot = time()
println("\n\n\t\tElapsed time tot: ", t_end_tot - t_start_tot, " seconds\n\n")
