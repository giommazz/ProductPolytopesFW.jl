# `get_solutions.jl`
using BPCGProduct
using FrankWolfe
using Printf

# filename = "intersecting_polytopes_n5_k10_mi1000000_v12-8-25-50-30-23-60-37-25-10_t20240516194427.jld2"
filename = "intersecting_polytopes_n3_k3_mi1000000_v12-8-25_t20240516194255.jld2"
n, k = extract_n_k_iters_from_filename(filename)

# Use parameters from YAML file
config = Config("examples/config.yml"; n=n, k=k, max_iterations=5*1000000) 

# Use ConvexHullOracle
cvxhflag = true

# Retrieve LMOs from previously generated instances
lmo_list_nonintersecting, lmo_list_intersecting = get_lmos(config, filename, cvxhflag=cvxhflag)

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
    
    # Gather data to plot
    push!(trajectories, trajectories_curr)
end

println("°°°°°°°°°°°°°°°°°°°°°°°°°°°")
println(length(trajectories))
println("°°°°°°°°°°°°°°°°°°°°°°°°°°°")

# Get, from each tuple: Iteration, Primal, Dual, Dual Gap, Time
function get_trajdata(iter::Tuple)
    return iter[1], iter[2], iter[3], iter[4], iter[5]
end

function print_trajdata(trajdata::Vector{Any}, print_iter::Int64, opt::Float64)
    # `trajdata[1][1][i]` = (Iteration, Primal, Dual, Dual Gap, Time) → tuple w/data about i-th iteration of FW
    trajdata = trajdata[1][1]
    
    println("It      Primal              Primal Gap          Dual Gap            Time")
    println("-------------------------------------------------------------------------")
    
    # Print iterations
    for (idx, iter) in enumerate(trajdata)
        i, prim, _, dgap, time = get_trajdata(iter)
        pgap = max(0, prim - opt)  # Ensure pgap is not negative
        if idx == 1 || idx == length(trajdata) || idx % print_iter == 0
            @Printf.printf("%-8d %-18.6e %-18.6e %-18.6e %-18.6e\n", i, prim, pgap, dgap, time)
        end
    end
end

# Labels for the plots
labels = ["Non-Intersecting", "Intersecting"]
traj_data = [trajectories[1][1], trajectories[2][1]]
# Plot trajectories
plot_trajectories(traj_data, labels, xscalelog=true)

#print_trajdata(trajectories_nonintersecting, 10, 2.540537e+02)
