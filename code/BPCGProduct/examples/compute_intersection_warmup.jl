# This script is used to quickly "warm up" the JIT compilation (==precompile) by running the code on a small instance
# So, you know, do not touch it. 
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
        println("\n\n\n ----------> Full Block-coordinate Away FW (ours)")
        _, _, _, _, td_full_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), AwayStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_full_bc_afw, trajectories_ni, trajectories_i, primal)
    end

    return trajectories_ni, trajectories_i
end


# ---------------------------------------------------------------------------------
# Use parameters from YAML file
config = Config("examples/config.yml")
print_config(config)
config_warmup = modify_config(config, k=2, n=15)
print_config(config_warmup)

# Generate instances 
println("********************************************************")
println("Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config_warmup)

# Optimal solution
primal = primal - 1     # Numerical reasons
basename = generate_filename(config)

# Labels for the plots
labels = ["F-BC-AFW"]

# execute main
println("\n\n********************************************************")
println("Running FW on the instances")
println("********************************************************")
trajectories_ni, trajectories_i = main(config_warmup, vertices, shifted_vertices, primal, labels, basename)