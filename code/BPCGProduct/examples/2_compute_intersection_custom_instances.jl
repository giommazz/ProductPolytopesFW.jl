# `compute_intersection_custom_instances.jl`
using BPCGProduct
using FrankWolfe

# Select instance file: contains two instances, with k polytopes, intersecting and non-intersecting
filename = "intersecting_polytopes_n10_k2_v20-21_t20240524113914.jld2"

n, k = extract_n_k_from_filename(filename)

# Use parameters from YAML file
config = Config("examples/config.yml"; n=n, k=k)
print_config(config)

# Use FrankWolfe.ConvexHullOracle LMOs (true) or FrankWolfe.MathOptLMO LMOs (false)
cvxhflag = false

# Load data and transform to Polyhedra.Polyhedron or JuMP.Model
vertices, shifted_vertices, primal, fw_gap = load_intersecting_polytopes(filename)    

# Retrieve nonintersecting and intersecting LMOs from previously generated instances
lmo_list = create_lmos(config, [vertices, shifted_vertices], cvxhflag=cvxhflag)

# run FW
trajectories = []

for (i, lmos) in enumerate(lmo_list)
    # nonintersecting flag
    ni_flag = i == 1
    println()
    if ni_flag println("\nNon intersecting") else println("\nIntersecting") end
    # Find possible subsets of size `config.k`
    # lmo_products = unique_combinations(config, lmo_list)
    
    prod_lmo = create_product_lmo(config, lmos)
    # println("\n\n\n---------------------------------------------------------")
    # println("LMOs in ProductLMO: ")
    # for i in [typeof(prod_lmo.lmos[i]) for i in 1:config.k] println("\t$i") end
    # println("---------------------------------------------------------")
    # Block-coordinate BPCG with CyclicUpdate
    println("----------> Cyclic Block-coordinate BPCG")
    _, _, _, _, trajectory_data = run_FW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
    
    # add optimal solution
    opt = 0.0
    if ni_flag opt = primal end

    # Replace "Primal" in the FW log with "Primal Gap", i.e. f(x) with f(x) - f(x*) 
    trajectory_data_pg = compute_primal_gap(trajectory_data, opt)
    
    # Gather data to plot
    append!(trajectories, trajectory_data_pg)
end

# Labels for the plots
labels = ["Non-Intersecting", "Intersecting"]
# Plot trajectories
plot_trajectories(trajectories, labels, yscalelog=false, xscalelog=true)

#print_trajdata(trajectories, 10, 2.540537e+02)