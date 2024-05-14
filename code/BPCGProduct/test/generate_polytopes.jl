# `generate_polytopes.jl`
using BPCGProduct
using BPCGProduct: Polyhedra
using BPCGProduct: JuMP

config = Config("test/config.yml") # Use parameters from YAML file

vertices, polytopes, intersecting_polytopes_polyhedra, intersecting_polytopes_jump, min_distance = generate_intersecting_polytopes(config)