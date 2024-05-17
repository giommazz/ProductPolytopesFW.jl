# `get_solutions.jl`
using BPCGProduct

config = Config("test/config.yml") # Use parameters from YAML file
filename = "intersecting_polytopes_n2_k10_v12-8-25-50-30-23-60-37-25-10_t20240516194119.jld2"
cvxflag = true

lmo_list_nonintersecting, lmo_list_intersecting = get_lmos(config, filename, cvxflag=cvxflag)
trajectories = get_solutions(lmo_list)