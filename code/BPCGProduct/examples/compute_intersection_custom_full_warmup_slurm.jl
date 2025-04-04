# `compute_intersection_custom_instances.jl`
# Run either within the Julia REPL as include("/examples/compute_intersection_custom_full.jl")
# or from linux terminal with: `julia --project=. examples/compute_intersection_custom_full.jl > test.log 2>&1`
#
# This includes a little script that warms up the REPL: it is executed before running the main script, so that when the main script is run, 
#   the package `BPCGProduct` has already been precompiled and the plots obtained by `main` don't show initial arbitrary overhead
using BPCGProduct
using FrankWolfe
using Plots





# ---------------------------------------------------------------------------------
# YAML PARAMETERS
config = Config("examples/config.yml")
print_config(config)
config_warmup = modify_config(config, k=2, n=15)
print_config(config_warmup)





# ---------------------------------------------------------------------------------
# MAIN FUNCTIONS
# Solve small instance (using `config_warmup`) to "warm-up" the REPL: this compiles `BPCGProduct`, so that no compilation needed upon running `main`
function repl_warmup(config::Config, vertices, shifted_vertices, primal, labels, basename)

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
        println("\n\n\n ----------> Full Block-coordinate Away FW (ours)")
        _, _, _, _, td_full_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), AwayStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_full_bc_afw, trajectories_ni, trajectories_i, primal)
    end

    return trajectories_ni, trajectories_i
end


function main(config::Config, path_to_results, vertices, shifted_vertices, primal, labels, basename)

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
        push_to_trajectories!(ni_flag, td_cyc_bc_fw, trajectories_ni, trajectories_i, primal)
        
        println("\n\n\n ----------> Cyclic Block-coordinate AFW")
        _, _, _, _, td_cyc_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), AwayStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_cyc_bc_afw, trajectories_ni, trajectories_i, primal)

        println("\n\n\n ----------> Stochastic Block-coordinate vanilla FW")
        _, _, _, _, td_stoc_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.StochasticUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_stoc_bc_fw, trajectories_ni, trajectories_i, primal)
        
        println("\n\n\n ----------> Stochastic Block-coordinate AFW")
        _, _, _, _, td_stoc_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.StochasticUpdate(), AwayStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_stoc_bc_afw, trajectories_ni, trajectories_i, primal)
        
        # println("\n\n\n ----------> Cyclic Block-coordinate BPFW")
        # _, _, _, _, td_cyc_bc_bpfw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_cyc_bc_bpfw, trajectories_ni, trajectories_i, primal)

        # ***************************************
        # Full block-coordinate methods
        println("\n\n\n ----------> Full Block-coordinate vanilla FW")
        _, _, _, _, td_full_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_full_bc_fw, trajectories_ni, trajectories_i, primal)

        println("\n\n\n ----------> Full Block-coordinate Away FW (ours)")
        _, _, _, _, td_full_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), AwayStep(), prod_lmo)  
        push_to_trajectories!(ni_flag, td_full_bc_afw, trajectories_ni, trajectories_i, primal)

        # println("\n\n\n ----------> Full Block-coordinate Blended Pairwise FW (ours)")
        # _, _, _, _, td_full_bc_bpfw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), FrankWolfe.BPCGStep(), prod_lmo)  
        # push_to_trajectories!(ni_flag, td_full_bc_bpfw, trajectories_ni, trajectories_i, primal)

        # ***************************************
        # Full methods
        println("\n\n\n ----------> Full FW")
        _, _, _, _, td_full_fw = run_FullFW(config, FrankWolfe.frank_wolfe, prod_lmo)    
        push_to_trajectories!(ni_flag, td_full_fw, trajectories_ni, trajectories_i, primal)

        println("\n\n\n ----------> Full AFW")
        _, _, _, _, td_full_afw = run_FullFW(config, FrankWolfe.away_frank_wolfe, prod_lmo)    
        push_to_trajectories!(ni_flag, td_full_afw, trajectories_ni, trajectories_i, primal)

        # println("\n\n\n ----------> Full BPFW")
        # _, _, _, _, td_full_bpfw = run_FullFW(config, FrankWolfe.blended_pairwise_conditional_gradient, prod_lmo)    
        # push_to_trajectories!(ni_flag, td_full_bpfw, trajectories_ni, trajectories_i, primal)
        
        # println("\n\n\n ----------> AP")
        # _, _, _, _, td_ap = run_AlternatingProjections(config, prod_lmo, true)    
        # # `FrankWolfe.alternating_projections` computes ||x-y|| rather than 1/2 ||x-y||
        # push_to_trajectories!(ni_flag, td_ap, trajectories_ni, trajectories_i, 2*primal)

        # Save trajectories
        # save_trajectories("examples/traj_$basename.jld2", trajectories_ni, trajectories_i)
    end

    log_data(trajectories_i, labels, path_to_results*"/times_i_"*basename)
    log_data(trajectories_ni, labels, path_to_results*"/times_ni_"*basename)

    return trajectories_ni, trajectories_i
end





# ---------------------------------------------------------------------------------
# WARM-UP SCRIPT
# Generate instances 
println("********************************************************")
println("Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config_warmup)
# Optimal solution
primal = primal - 1     # Numerical reasons
basename = generate_filename(config_warmup)
# Labels for the plots
labels = ["F-BC-AFW"]
# execute main
println("\n\n********************************************************")
println("Running FW on the instances")
println("********************************************************")
trajectories_ni, trajectories_i = repl_warmup(config_warmup, vertices, shifted_vertices, primal, labels, basename)





# ---------------------------------------------------------------------------------
# MAIN SCRIPT
results_directory = "examples/results_linesearch_afw" # "results_shortstep", 

# Generate instances 
println("********************************************************")
println("Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config)
# Optimal solution
primal = primal - 1     # Numerical reasons
basename = generate_filename(config)

# Labels for the plots
labels = ["C-BC-FW", "C-BC-AFW", "S-BC-FW", "S-BC-AFW", "F-BC-FW", "F-BC-AFW", "F-FW", "F-AFW"] # ["C-BC-FW", "C-BC-AFW", "C-BC-BPFW", "F-BC-FW", "F-BC-AFW", "F-BC-BPFW", "F-FW", "F-AFW", "F-BPFW", "AP"]

# execute main
println("\n\n********************************************************")
println("Running FW on the instances")
println("********************************************************")
trajectories_ni, trajectories_i = main(config, results_directory, vertices, shifted_vertices, primal, labels, basename)

# Plot trajectories
fig_ni_filename = results_directory*"/plot_ni_$basename"
fig_i_filename = results_directory*"/plot_i_$basename"
# Generate plots but do not pass `filename` argument (so .png is not automatically saved)
fig_ni = plot_trajectories(trajectories_ni, labels, yscalelog=true, xscalelog=true)
fig_i = plot_trajectories(trajectories_i, labels, yscalelog=true, xscalelog=true)

# Manually save only the PDF versions
#Plots.plot!(fig_ni, size=(1200, 800))  # Larger figure size
#Plots.plot!(fig_i, size=(1200, 800))  # Larger figure size
Plots.savefig(fig_ni, fig_ni_filename*".pdf")  
Plots.savefig(fig_i, fig_i_filename*".pdf")