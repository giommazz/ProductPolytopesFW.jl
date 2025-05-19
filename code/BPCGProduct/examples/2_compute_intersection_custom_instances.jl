# `2_compute_intersection_custom_instances.jl`
using ProductPolytopesAFW
using FrankWolfe

function main(config::Config)

    # Load data and transform to Polyhedra.Polyhedron or JuMP.Model
    vertices, shifted_vertices, primal, _ = load_polytopes(filename)    

    # Retrieve nonintersecting and intersecting LMOs from previously generated instances
    lmo_list = create_lmos(config, [vertices, shifted_vertices])

    # Will contain data about diafferent FW runs, for non-intersecting and intersecting polytopes
    trajectories_ni, trajectories_i = [], []

    for (i, lmos) in enumerate(lmo_list)
        # nonintersecting flag
        ni_flag = i == 1
        println()
        if ni_flag println("\nNon intersecting") else println("\nIntersecting") end

        prod_lmo = create_product_lmo(config, lmos)

        # Run Frank-Wolfe algorithms and alternating projectiopng, then record trajectory data
        println("\n\n\n ----------> Cyclic Block-coordinate vanilla CG")
        _, _, _, _, td_cyc_bc_cg = run_FW(config, FrankWolfe.CyclicUpdate(), prod_lmo)
        println("\n\n\n ----------> Cyclic Block-coordinate BPCG")
        _, _, _, _, td_cyc_bc_bpcg = run_FW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
        println("\n\n\n ----------> Full Block-coordinate BPCG")
        _, _, _, _, td_full_bc_cg = run_FW(config, FrankWolfe.FullUpdate(), FrankWolfe.BPCGStep(), prod_lmo)  
        println("\n\n\n ----------> Full BPCG")
        _, _, _, _, td_bpcg = run_FW(config, prod_lmo)    
        println("\n\n\n ----------> AP")
        _, _, _, _, td_ap = run_FW(config, prod_lmo, true)
        
        push_to_trajectories!(ni_flag, td_cyc_bc_cg, trajectories_ni, trajectories_i, primal)
        push_to_trajectories!(ni_flag, td_cyc_bc_bpcg, trajectories_ni, trajectories_i, primal)
        push_to_trajectories!(ni_flag, td_full_bc_cg, trajectories_ni, trajectories_i, primal)
        push_to_trajectories!(ni_flag, td_bpcg, trajectories_ni, trajectories_i, primal)
        push_to_trajectories!(ni_flag, td_ap, trajectories_ni, trajectories_i, primal)

        # Save trajectories
        save_trajectories("examples/traj_$basename.jld2", trajectories_ni, trajectories_i)

    end

    return trajectories_ni, trajectories_i
end


# ---------------------------------------------------------------------------------

# Select instance file: contains two instances, with k polytopes, intersecting and non-intersecting
#filename = "intersecting_polytopes_n5_k2_v10-11_t20240524121720.jld2" 
#filename = "intersecting_polytopes_n10_k2_v20-21_t20240524152154.jld2" 
filename = "intersecting_polytopes_n50_k2_v94-95_t20240524131142.jld2"
basename = base_name(filename)

n, k = extract_n_k_from_filename(filename)

# Use parameters from YAML file
config = Config("examples/config.yml"; n=n, k=k)
print_config(config)

# execute main
trajectories_ni, trajectories_i = main(config)

# Labels for the plots
labels = ["C-BC-FW", "C-BC-BPCG", "F-BC-BPCG", "F-BPCG", "AP"]
# Plot trajectories
plot_trajectories(trajectories_ni, labels, yscalelog=false, xscalelog=true, filename="examples/plot_ni_$basename.png")
plot_trajectories(trajectories_i, labels, yscalelog=false, xscalelog=true, filename="examples/plot_i_$basename.png")
