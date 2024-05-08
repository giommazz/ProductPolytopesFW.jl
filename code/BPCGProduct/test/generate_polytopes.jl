# `generate_polytopes.jl`
using BPCGProduct

config = Config("test/config.yml") # Use parameters from YAML file

function main(config)
    model_A, x_A, _, _ = generate_polytope(2, config)
    println(model_A)
    direction = [rand(-1.0:1.0) for _ in 1:config.n]
    
    # Define a direction to find a vertex
    println(direction)
    
    # Find a vertex inside polytope A
    vertex_A = find_vertex_in_polytope(model_A, x_A, direction)
    vertex_B = find_vertex_in_polytope(model_A, x_A, -direction)
    println("Vertex in Polytope A:", vertex_A)
    println("Opposite Vertex in Polytope A:", vertex_B)

    model_B, translated_x_B = setup_translated_polytope_B(config, vertex_A)
    println(model_B)
    println(translated_x_B)

    #poly, fig = polytope_from_jump(m1)
end

main(config)