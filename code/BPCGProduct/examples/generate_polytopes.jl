# `generate_polytopes.jl`
using BPCGProduct

# Use parameters from YAML file
config = Config("examples/config.yml")
# Generate instances 
vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config)

# Save data to .jld file
filename = save_intersecting_polytopes(config, vertices, shifted_vertices, primal, fw_gap)
# As an alternative, one can give a filename
# filename = "example.jdl2"
# save_intersecting_polytopes(filename, vertices, shifted_vertices, primal, fw_gap)