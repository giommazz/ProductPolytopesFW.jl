# `get_solutions.jl`
using BPCGProduct

config = Config("examples/config.yml") # Use parameters from YAML file
filename = "intersecting_polytopes_n2_k10_v12-8-25-50-30-23-60-37-25-10_t20240516194119.jld2"
# Use ConvexHullOracle
cvxflag = false

# Build LMOs from previously generated instances
lmo_list_nonintersecting, lmo_list_intersecting = get_lmos(config, filename, cvxflag=cvxflag)

# run FW
trajectories_nonintersecting = get_solutions(config, lmo_list_nonintersecting)
trajectories_intersecting = get_solutions(config, lmo_list_intersecting)



#TODO: AUMENTA IL NUMERO DI ITERAZIONI!!!!