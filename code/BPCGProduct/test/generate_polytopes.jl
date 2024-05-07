# `generate_polytopes.jl`
using BPCGProduct
using BPCGProduct: CDDLib
using BPCGProduct: Polyhedra
using BPCGProduct: Makie


config = Config("test/config.yml") # Use parameters from YAML file

function main(config)
    model_A, _ = generate_polytope(4, config)
    poly = Polyhedra.polyhedron(model_A, CDDLib.Library(:exact))
    Makie.mesh(Polyhedra.Mesh(poly), color=:blue)
end

main(config)