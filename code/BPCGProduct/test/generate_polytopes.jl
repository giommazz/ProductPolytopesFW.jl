# `generate_polytopes.jl`
using BPCGProduct

# Use parameters from YAML file
config = Config("test/config.yml")
# Generate instances 
vertices, shifted_vertices, polytopes, intersecting_polytopes_polyhedra, intersecting_polytopes_jump = 
    generate_intersecting_polytopes(config)

# Save data to .jld file
filename = save_intersecting_polytopes(config, vertices, shifted_vertices, intersecting_polytopes_jump)
# As an alternative, one can give a filename
# filename = "prova.jdl2"
# save_intersecting_polytopes(filename, vertices, polytopes, intersecting_polytopes_polyhedra, intersecting_polytopes_jump)

# Load the intersecting polytopes data for further analysis
#verticesj, polytopesj, intersecting_polytopes_polyhedraj, intersecting_polytopes_jumpj = 
# TODO: AGGIUSTA
verticesj = load_intersecting_polytopes(filename)

# println(typeof(verticesj))
# println(typeof(polytopesj))
# println(typeof(intersecting_polytopes_polyhedraj))
# println(typeof(intersecting_polytopes_jumpj))