# `get_solutions.jl`
using BPCGProduct
using FrankWolfe

# filename = "intersecting_polytopes_n5_k10_mi1000000_v12-8-25-50-30-23-60-37-25-10_t20240516194427.jld2"
filename = "intersecting_polytopes_n3_k3_mi1000000_v12-8-25_t20240516194255.jld2"
n, k = extract_n_k_iters_from_filename(filename)

# Use parameters from YAML file
config = Config("examples/config.yml"; n=n, k=k, max_iterations=5*1000000) 

# Use ConvexHullOracle
cvxhflag = false

# Retrieve LMOs from previously generated instances
lmo_list_nonintersecting, lmo_list_intersecting = get_lmos(config, filename, cvxhflag=cvxhflag)

# run FW
trajectories_nonintersecting, trajectories_intersecting = [], []

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
    if ni_flag
        push!(trajectories_nonintersecting, trajectories_curr)
    else
        push!(trajectories_intersecting, trajectories_curr)
    end
end

println("°°°°°°°°°°°°°°°°°°°°°°°°°°°")
println(typeof(trajectories_intersecting))
println("°°°°°°°°°°°°°°°°°°°°°°°°°°°")
println(trajectories_intersecting)
println("°°°°°°°°°°°°°°°°°°°°°°°°°°°")
println(trajectories_intersecting[1])
println("°°°°°°°°°°°°°°°°°°°°°°°°°°°")
println(trajectories_intersecting[1][1])
println("°°°°°°°°°°°°°°°°°°°°°°°°°°°")
println(trajectories_intersecting[1][1][1])

plot_trajectories(trajectories_intersecting[1], "I")
#plot_trajectories([trajectories_intersecting, trajectories_nonintersecting], ["I", "NI"])#, xscalelog=true)