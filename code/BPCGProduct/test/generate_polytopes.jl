# `generate_polytopes.jl`
using BPCGProduct

config = Config("test/config.yml") # Use parameters from YAML file

function main(config)
    m1, _ = generate_polytope(2, config)
    println(m1)
    m2, _ = generate_polytope(config)
    println()
    println(m2)
    poly, fig = polytope_from_jump(m1)
end

main(config)