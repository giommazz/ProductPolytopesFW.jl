# `compute_intersection_custom_instances.jl`
using BPCGProduct
using FrankWolfe

function main(config::Config, vertices, shifted_vertices, primal, labels, basename)

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
        
        println("\n\n\n ----------> Cyclic Block-coordinate vanilla FW")
        _, _, _, _, td_cyc_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_cyc_bc_fw, trajectories_ni, trajectories_i, primal)

        # println("\n\n\n ----------> Cyclic Block-coordinate BPFW")
        # _, _, _, _, td_cyc_bc_bpcd = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_cyc_bc_bpcg, trajectories_ni, trajectories_i, primal)

        println("\n\n\n ----------> Full Block-coordinate Blended Pairwise FW (ours)")
        _, _, _, _, td_full_bc_bpcg = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), FrankWolfe.BPCGStep(), prod_lmo)  
        push_to_trajectories!(ni_flag, td_full_bc_bpcg, trajectories_ni, trajectories_i, primal)

        println("\n\n\n ----------> Full Block-coordinate Away FW (ours)")
        _, _, _, _, td_full_bc_acg = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), AwayStep(), prod_lmo)  
        push_to_trajectories!(ni_flag, td_full_bc_acg, trajectories_ni, trajectories_i, primal)

        println("\n\n\n ----------> Full BPFW")
        _, _, _, _, td_full_bpcg = run_FullBlendedPairwiseFW(config, prod_lmo)    
        push_to_trajectories!(ni_flag, td_full_bpcg, trajectories_ni, trajectories_i, primal)
        
        # println("\n\n\n ----------> AP")
        # _, _, _, _, td_ap = run_AlternatingProjections(config, prod_lmo, true)    
        # # `FrankWolfe.alternating_projections` computes ||x-y|| rather than 1/2 ||x-y||
        # push_to_trajectories!(ni_flag, td_ap, trajectories_ni, trajectories_i, 2*primal)

        # Save trajectories
        # save_trajectories("examples/traj_$basename.jld2", trajectories_ni, trajectories_i)
    end

    log_data(trajectories_i, labels, "examples/"*results_directory*"/times_i_"*basename)
    log_data(trajectories_ni, labels, "examples/"*results_directory*"/times_ni_"*basename)

    return trajectories_ni, trajectories_i
end


# ---------------------------------------------------------------------------------
# Use parameters from YAML file
config = Config("examples/config.yml")
print_config(config)

results_directory = "results_linesearch_acg" # "results_shortstep", 

# Generate instances 
println("********************************************************")
println("Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config)

# Optimal solution
primal = primal - 1     # Numerical reasons
basename = generate_filename(config)

# Labels for the plots
labels = ["C-BC-FW", "F-BC-BPFW", "F-BC-AFW", "F-BPFW"] # ["C-BC-FW", "C-BC-BPFW", "F-BC-BPFW", "F-BC-AFW", "F-BPFW", "AP"]

# execute main
println("\n\n********************************************************")
println("Running FW on the instances")
println("********************************************************")
trajectories_ni, trajectories_i = main(config,vertices, shifted_vertices, primal, labels, basename)


# Plot trajectories
plot_trajectories(trajectories_ni, labels, yscalelog=true, xscalelog=true, filename="examples/"*results_directory*"/plot_ni_$basename.png")
plot_trajectories(trajectories_i, labels, yscalelog=true, xscalelog=true, filename="examples/"*results_directory*"/plot_i_$basename.png")