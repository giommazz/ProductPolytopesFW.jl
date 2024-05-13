# `generate_polytopes.jl`
using BPCGProduct
using BPCGProduct: Polyhedra
using BPCGProduct: JuMP

config = Config("test/config.yml") # Use parameters from YAML file

generate_intersecting_polytopes(config)