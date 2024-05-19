# `get_solutions.jl`
using BPCGProduct

filename = "intersecting_polytopes_n2_k10_v12-8-25-50-30-23-60-37-25-10_t20240516194119.jld2"
n, k = extract_n_k_from_filename(filename)
# Use parameters from YAML file
config = Config("examples/config.yml", Dict("n"=> n, "k" => k)) 
println(config)
readline()

# Use ConvexHullOracle
cvxflag = false
config = update_config_with_n_k(config, filename)
# Build LMOs from previously generated instances
lmo_list_nonintersecting, lmo_list_intersecting = get_lmos(config, filename, cvxflag=cvxflag)

# run FW
trajectories_nonintersecting = get_solutions(config, lmo_list_nonintersecting)
trajectories_intersecting = get_solutions(config, lmo_list_intersecting)



#TODO: AUMENTA IL NUMERO DI ITERAZIONI!!!!