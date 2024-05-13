# `generate_polytopes.jl`
using BPCGProduct
using BPCGProduct: Polyhedra
using BPCGProduct: JuMP

config = Config("test/config.yml") # Use parameters from YAML file

# Main function to generate and move polytopes
function main(config)
    # Number of points for each polytope
    n_points = [4, 7]

    # Generate random, non-intersecting bounds
    bounds_list = generate_non_intersecting_bounds(config)

    # Initialize empty lists for vertices, Polyhedra polytopes, and Polyhedra intersecting polytopes, JuMP intersecting polytopes
    vertices = Vector{Matrix{Float64}}()
    polytopes = Vector{Polyhedra.Polyhedron}()
    intersecting_polytopes_polyhedra = Vector{Polyhedra.Polyhedron}()
    intersecting_polytopes_jump = Vector{JuMP.Model}()

    # Generate k polytopes within the specified bounds
    for i in 1:config.k
        # Generate vertices and corresponding polytope
        verts, poly = generate_polytope(config, n_points[i], bounds_list[i])
        
        # Append to `vertices` and `polytopes`
        push!(vertices, verts)
        push!(polytopes, poly)
    end

    # Add first polytope to `intersecting_polytopes_polyhedra` and to `intersecting_polytopes_jump`
    push!(intersecting_polytopes_polyhedra, polytopes[1])
    push!(intersecting_polytopes_jump, polyhedra_to_jump(config, polytopes[1]))

    # Move each subsequent polytope to intersect with the first one
    for i in 2:config.k
        # Move iᵗʰ polytope so that it intersects with first one
        intersecting_polytope = intersect_polytopes(config, vertices[1], vertices[i], polytopes[i])
        # Append to the list of intersecting polytopes
        push!(intersecting_polytopes_polyhedra, intersecting_polytope)
        push!(intersecting_polytopes_jump, polyhedra_to_jump(config, intersecting_polytope))
    end

    for i in 1:config.k
        println("°°°°°°°°°°°°°°°°---------------------------------")
        println(intersecting_polytopes_polyhedra[i])
        println(intersecting_polytopes_jump[i])
        for (_, constraint) in JuMP.all_constraints(intersecting_polytopes_jump[i])
            println("Constraint _: $constraint")
        end
        println()
    end



    println(Polyhedra.npoints(intersect(intersecting_polytopes_polyhedra[1], intersecting_polytopes_polyhedra[2])))

    

    #plot_polytopes(intersecting_polytopes_polyhedra)
end

main(config)