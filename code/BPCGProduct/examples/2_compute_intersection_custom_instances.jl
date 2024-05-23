# `compute_intersection_custom_instances.jl`
using BPCGProduct
using FrankWolfe

# Select instance file: contains two instances, with k polytopes, intersecting and non-intersecting
filename = "intersecting_polytopes_n10_k2_v20-21_t20240523172111.jld2"

n, k = extract_n_k_from_filename(filename)
max_iterations = 5*1000000

# Use parameters from YAML file
config = Config("examples/config.yml"; n=n, k=k, max_iterations=max_iterations)

# Use FrankWolfe.ConvexHullOracle LMOs (true) or FrankWolfe.MathOptLMO LMOs (false)
cvxhflag = false

# Load data and transform to Polyhedra.Polyhedron or JuMP.Model
vertices, shifted_vertices, primal, fw_gap = load_intersecting_polytopes(filename)    

# Retrieve LMOs from previously generated instances
lmo_list_nonintersecting, lmo_list_intersecting = create_lmos(config, vertices, shifted_vertices, cvxhflag=cvxhflag)

# run FW
trajectories = []

for lmo_list in [lmo_list_nonintersecting, lmo_list_intersecting]
    ni_flag = lmo_list == lmo_list_nonintersecting
    println()
    if ni_flag println("\nNon intersecting") else println("\nIntersecting") end
    # Find possible subsets of size `config.k`
    # lmo_products = unique_combinations(config, lmo_list)
    
    prod_lmo = create_product_lmo(config, lmo_list)
    # println("\n\n\n---------------------------------------------------------")
    # println("LMOs in ProductLMO: ")
    # for i in [typeof(prod_lmo.lmos[i]) for i in 1:config.k] println("\t$i") end
    # println("---------------------------------------------------------")
    # Block-coordinate BPCG with CyclicUpdate
    println("----------> Cyclic Block-coordinate BPCG")
    trajectories_curr = run_FW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
    
    # add optimal solution
    opt = 0.0
    if ni_flag opt = primal end

    # Replace "Primal" in the FW log with "Primal Gap", i.e. f(x) with f(x) - f(x*) 
    trajectories_curr_pg = compute_primal_gap(trajectories_curr, opt)
    
    # Gather data to plot
    append!(trajectories, trajectories_curr_pg)
end

# Labels for the plots
labels = ["Non-Intersecting", "Intersecting"]
# Plot trajectories
plot_trajectories(trajectories, labels, yscalelog=false, xscalelog=true)

#print_trajdata(trajectories, 10, 2.540537e+02)