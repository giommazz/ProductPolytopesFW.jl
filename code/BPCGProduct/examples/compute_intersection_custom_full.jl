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
        
        println("\n\n\n ----------> Cyclic Block-coordinate vanilla CG")
        _, _, _, _, td_cyc_bc_fw = run_FW(config, FrankWolfe.CyclicUpdate(), prod_lmo)
        push_to_trajectories!(ni_flag, td_cyc_bc_fw, trajectories_ni, trajectories_i, primal)

        # println("\n\n\n ----------> Cyclic Block-coordinate BPCG")
        # _, _, _, _, td_cyc_bc_bpcg = run_FW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_cyc_bc_bpcg, trajectories_ni, trajectories_i, primal)

        println("\n\n\n ----------> Full Block-coordinate BPCG")
        _, _, _, _, td_full_bc_cg = run_FW(config, FrankWolfe.FullUpdate(), FrankWolfe.BPCGStep(), prod_lmo)  
        push_to_trajectories!(ni_flag, td_full_bc_cg, trajectories_ni, trajectories_i, primal)

        # println("\n\n\n ----------> Full BPCG")
        # _, _, _, _, td_bpcg = run_FW(config, prod_lmo)    
        # push_to_trajectories!(ni_flag, td_bpcg, trajectories_ni, trajectories_i, primal)
        
        println("\n\n\n ----------> AP")
        _, _, _, _, td_ap = run_FW(config, prod_lmo, true)    
        push_to_trajectories!(ni_flag, td_ap, trajectories_ni, trajectories_i, primal)

        # Save trajectories
        # save_trajectories("examples/traj_$basename.jld2", trajectories_ni, trajectories_i)

    end
    
    log_data(trajectories_i, labels, "times_i_"*basename)
    log_data(trajectories_ni, labels, "times_ni_"*basename)

    return trajectories_ni, trajectories_i
end


# ---------------------------------------------------------------------------------
# Use parameters from YAML file
config = Config("examples/config.yml")
print_config(config)

# Generate instances 
println("********************************************************")
println("Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config)
# Numerical reasons
primal = primal - 1
basename = generate_filename(config, vertices)

# Labels for the plots
labels = ["C-BC-FW", "F-BC-BPCG", "AP"]# ["C-BC-FW", "C-BC-BPCG", "F-BC-BPCG", "F-BPCG", "AP"]

# execute main
println("\n\n********************************************************")
println("Running FW on the instances")
println("********************************************************")
trajectories_ni, trajectories_i = main(config,vertices, shifted_vertices, primal, labels, basename)


# Plot trajectories
plot_trajectories(trajectories_ni, labels, yscalelog=true, xscalelog=true, filename="plot_ni_$basename.png")
plot_trajectories(trajectories_i, labels, yscalelog=true, xscalelog=true, filename="plot_i_$basename.png")
