# `get_solutions.jl`
using BPCGProduct
using FrankWolfe

filename = "intersecting_polytopes_n3_k3_v12-8-25_t20240516194255.jld2"
n, k = extract_n_k_from_filename(filename)
# Use parameters from YAML file
config = Config("examples/config.yml"; n=n, k=k) 

# Use ConvexHullOracle
cvxhflag = false
config = update_config_with_n_k(config, filename)

# Retrieve LMOs from previously generated instances
lmo_list_nonintersecting, lmo_list_intersecting = get_lmos(config, filename, cvxhflag=cvxhflag)

# run FW
trajectories_nonintersecting, trajectories_intersecting = [], []

for lmo_list in [lmo_list_nonintersecting, lmo_list_intersecting]
    # Find possible subsets of size `config.k`
    lmo_products = unique_combinations(config, lmo_list)
    for lmos in lmo_products
        prod_lmo = create_product_lmo(config, lmos)
        println("\n\n\n---------------------------------------------------------")
        println("LMOs: ", [typeof(prod_lmo.lmos[i]) for i in 1:config.k])
        println("---------------------------------------------------------")
        # Block-coordinate BPCG with CyclicUpdate
        println("\n\n\n ----------> Cyclic Block-coordinate BPCG")
        trajectories_curr = run_FW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
        if lmo_list == lmo_list_nonintersecting
            push!(trajectories_nonintersecting, trajectories_curr)
        else
            push!(trajectories_intersecting, trajectories_curr)
        end
    end
end





# TODO: AUMENTA IL NUMERO DI ITERAZIONI!!!!