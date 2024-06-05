# `generate_polytopes.jl`
using BPCGProduct

# Use parameters from YAML file
config = Config("examples/config.yml")
print_config(config)

# Generate instances 
t = @elapsed vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config)
println("°°°°°°°°°°°°°°°°°°°°°°°°°°°° time generate_polytopes: $t")

# Save data to .jld file
filename = save_polytopes(config, vertices, shifted_vertices, primal, fw_gap)
# As an alternative, one can give a filename
# filename = "example.jdl2"
# save_polytopes(filename, vertices, shifted_vertices, primal, fw_gap)