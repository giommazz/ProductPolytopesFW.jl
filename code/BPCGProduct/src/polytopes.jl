# `polytopes.jl`

# Function to generate non-overlapping bounds for multiple dimensions
function generate_non_intersecting_bounds(config::Config; margin::Float64 = 1.0)
    bounds_list = Vector{Vector{Tuple{Float64, Float64}}}(undef, config.k)

    # Generate the first bounds randomly
    bounds_list[1] = [(rand(0.0:5.0), rand(5.0:10.0)) for _ in 1:config.n]

    # Generate subsequent bounds while ensuring no overlap
    for i in 2:config.k
        bounds = Vector{Tuple{Float64, Float64}}(undef, config.n)

        for d in 1:config.n
            prev_min, prev_max = bounds_list[i - 1][d]
            lower_bound = prev_max + margin
            upper_bound = lower_bound + rand(5.0:10.0)
            bounds[d] = (lower_bound, upper_bound)
        end

        bounds_list[i] = bounds
    end

    return bounds_list
end

# Function to generate a random matrix with bounds for each column
function generate_polytope(config::Config, points::Int, bounds::Vector{Tuple{Float64, Float64}})
    # Check that the number of columns specified in `config.n` matches the length of `bounds`
    if config.n != length(bounds)
        error("Mismatch: The number of columns in `config` (config.n = $(config.n)) must match the number of elements in `bounds` (length = $(length(bounds))).")
    end

    # Initialize a matrix of zeros with the given number of points (rows) and columns
    vertices = zeros(Float64, points, config.n)

    # Fill each column based on the specified bounds
    for j in 1:config.n
        lower_bound, upper_bound = bounds[j]
        # Generate random numbers within the given bounds
        vertices[:, j] = lower_bound .+ (upper_bound - lower_bound) .* rand(Float64, points)
    end
    # Generate polytope as convex hull of given vertices
    poly = polyhedron(vrep(vertices), CDDLib.Library())
    return vertices, poly
end

# Function to find the closest points between two sets of points
function closest_pair(config::Config, V1::Matrix, V2::Matrix)

    # Check that both matrices have the correct number of columns
    if size(V1, 2) != config.n || size(V2, 2) != config.n
        error("Both point sets must have the dimension specified in `config.n`.")
    end

    # Number of points in each set
    n_V1 = size(V1, 1)
    n_V2 = size(V2, 1)

    # Initialize minimum distance to infinity and empty indices
    min_dist = Inf
    closest_indices = (-1, -1)

    # Iterate through all pairs and find the minimum distance
    for i in 1:n_V1
        for j in 1:n_V2
            dist = sum((V1[i, :] - V2[j, :]).^2)  # Squared Euclidean distance
            if dist < min_dist
                min_dist = dist
                closest_indices = (i, j)
            end
        end
    end

    # Retrieve the closest points
    v1_closest = V1[closest_indices[1], :]
    v2_closest = V2[closest_indices[2], :]

    return v1_closest, v2_closest, min_dist
end
# (Multiple Dispatch) Here the first input is a vector (single point)
function closest_pair(config::Config, v::Vector, V2::Matrix)
    # Check that the vector and matrix have the correct number of columns
    if length(v) != config.n || size(V2, 2) != config.n
        error("Both the point and the point set must have the dimension specified in `config.n`.")
    end

    # Number of points in V2
    n_V2 = size(V2, 1)

    # Initialize minimum distance to infinity and the index of the closest point
    min_dist = Inf
    closest_index = -1

    # Iterate through all points in V2 to find the closest one
    for j in 1:n_V2
        dist = sum((v - V2[j, :]).^2)  # Squared Euclidean distance
        if dist < min_dist
            min_dist = dist
            closest_index = j
        end
    end

    # Retrieve the closest point
    v2_closest = V2[closest_index, :]

    return v2_closest, min_dist
end

# Function to intersect two polytopes, moving the second onto the closest vertex of the first
function intersect_polytopes(config::Config, V1::Matrix{T}, V2::Matrix{T}, polytope_to_move::Polyhedron{T}) where T
    
    # Find the closest pair of points between the two sets of points
    v1_closest, v2_closest, distance = closest_pair(config, V1, V2)

    # Create a translation vector (offset) from p_closest - q_closest
    direction = v1_closest - v2_closest

    # Create "fake" offset polytope
    P = polyhedron(vrep([direction]), CDDLib.Library())
    # Shift `polytope_to_move` along the offset
    Pint = polytope_to_move + P

    # Find vertices of `Pint` (`Polyhedron.ext` gives the V-representation of a polytope)
    generators = Polyhedra.removevredundancy(Pint.ext, GLPK.Optimizer)
    # Remove points that are not vertices of `Pint`
    Polyhedra.removevredundancy!(Pint)

    # Return translated polytope and other info
    return Pint, generators, v1_closest, v2_closest, distance
end
# (Multiple Dispatch) Function to intersect two polytopes, moving the second unto a given vertex of the first
function intersect_polytopes(config::Config, v::Vector, V2::Matrix{T}, polytope_to_move::Polyhedron{T}) where T
    
    # Find the closest pair of points between the two sets of points
    v2_closest, distance = closest_pair(config, v, V2)

    # Create a translation vector (offset) from p_closest - q_closest
    direction = v - v2_closest

    # Create "fake" offset polytope
    P = polyhedron(vrep([direction]), CDDLib.Library())
    # Shift `polytope_to_move` along the offset
    Pint = polytope_to_move + P
    
    # Find vertices of `Pint` (`Polyhedron.ext` gives the V-representation of a polytope)
    generators = Polyhedra.removevredundancy(Pint.ext, GLPK.Optimizer)
    # Remove points that are not vertices of `Pint`
    Polyhedra.removevredundancy!(Pint)

    # Return translated polytope and other info
    return Pint, generators, v2_closest, distance
end

# Function to generate a polytope with a given JuMP model
function polyhedra_to_jump(config::Config, polytope::Polyhedra.Polyhedron{T}) where T
    model = Model()
    @variable(model, x[1:config.n])
    @constraint(model, x in polytope)
    
    return model
end

# Main function to generate and move polytopes
# `n_points` contains n. of vertices used to generate each polytope
function generate_intersecting_polytopes(config)

    println("Generating $(config.k) intersecting polytopes")
    
    # Generate random, non-intersecting bounds
    bounds_list = generate_non_intersecting_bounds(config)

    # Initialize empty lists for vertices, Polyhedra polytopes, and Polyhedra intersecting polytopes, JuMP intersecting polytopes
    vertices = Vector{Matrix{Float64}}()
    polytopes = Vector{Polyhedra.Polyhedron}()
    intersecting_polytopes_polyhedra = Vector{Polyhedra.Polyhedron}()
    intersecting_polytopes_jump = Vector{JuMP.Model}()
    min_distance = Float64

    # Generate k polytopes within the specified bounds
    for i in 1:config.k
        # Generate vertices and corresponding polytope
        verts, poly = generate_polytope(config, config.n_points[i], bounds_list[i])
        
        # Append to `vertices` and `polytopes`
        push!(vertices, verts)
        push!(polytopes, poly)
    end

    # Push P₁ to `intersecting_polytopes_polyhedra` and to `intersecting_polytopes_jump`
    push!(intersecting_polytopes_polyhedra, polytopes[1])
    push!(intersecting_polytopes_jump, polyhedra_to_jump(config, polytopes[1]))
    
    # Move P₂ towards P₁ so that they intersect (at least) in `v₁`
    intersecting_polytope, _, v₁, _, distance = intersect_polytopes(config, vertices[1], vertices[2], polytopes[2])
    push!(intersecting_polytopes_polyhedra, intersecting_polytope)
    push!(intersecting_polytopes_jump, polyhedra_to_jump(config, intersecting_polytope))
    intersection = Polyhedra.npoints(intersect(polytopes[1], intersecting_polytope))
    println("Intersection of P₁ and P₂ contains $intersection points")
    @assert intersection ≥ 1 "There must be at least one point in P₁ ∩ Pᵢ"
    
    # Move each subsequent polytope Pₖ, for k > 3, to intersect with P₁
    for i in 3:config.k
        # Move iᵗʰ polytope so that it intersects with the others (at least) in `v₁`
        intersecting_polytope, _, _, distance = intersect_polytopes(config, v₁, vertices[i], polytopes[i])

        # Append to list of intersecting polytopes
        push!(intersecting_polytopes_polyhedra, intersecting_polytope)
        push!(intersecting_polytopes_jump, polyhedra_to_jump(config, intersecting_polytope))
        intersection = Polyhedra.npoints(intersect(polytopes[1], intersecting_polytope))
        println("Intersection of P₁ and Pᵢ [for i=$i] contains $intersection points")
        @assert intersection ≥ 1 "There must be at least one point in P₁ ∩ Pᵢ"
    end

    return vertices, polytopes, intersecting_polytopes_polyhedra, intersecting_polytopes_jump, min_distance
end

# Compute 1/2 ∑ᵢ₌₁ᵏ⁻¹∑ⱼ₌ᵢ₊₁ᵏ || xⁱ - xʲ ||₂²
function compute_min_distance(config::Config, polytopes::Vector{Polyhedra.Polyhedron})
    
    min_dist = 0.0
    
    for i in 1:(config.k - 1)
        for j in (i + 1):config.k
            # TODO: USE min_distance_between_polytopes
            min_dist += dist
        end
    end
    
    return 0.5 * min_dist_sum
end

# Function to formulate and solve the minimum distance problem
function min_distance_between_polytopes(config::Config, poly1::Polyhedra.Polyhedron, poly2::Polyhedra.Polyhedron)
    
    model = Model(GLPK.Optimizer)
    
    @variable(model, p[1:config.n])
    @variable(model, q[1:config.n])
    
    @constraint(model, p in poly1)
    @constraint(model, q in poly2)
    
    @objective(model, Min, (1/2)*sum((p[i] - q[i])^2 for i in 1:config.n))
    
    optimize!(model)
    
    return objective_value(model)
end

# for i in 1:config.k
#     println("°°°°°°°°°°°°°°°°---------------------------------> [$i]")
    
#     println("\tVREP")
#     for constraint in JuMP.list_of_constraint_types(intersecting_polytopes_jump[i])
#         println("\t\tConstraint: $constraint")
#     end
#     println()
#     println("\t\tall constraints: ", JuMP.all_constraints(intersecting_polytopes_jump[i]; include_variable_in_set_constraints = true))
#     println("number of constraints: $(JuMP.num_constraints(intersecting_polytopes_jump[i]; count_variable_in_set_constraints = true))")
#     println()   
#     println()


#     PP = Polyhedra.hrep(intersecting_polytopes_polyhedra[i])
#     println("\tHREP")
#     PPJ = polyhedra_to_jump(config, intersecting_polytopes_polyhedra[i])
#     for constraint in JuMP.list_of_constraint_types(PPJ)
#         println("\t\tConstraint: $constraint")
#     end
#     println()
#     println("\t\tall constraints: ", JuMP.all_constraints(intersecting_polytopes_jump[i]; include_variable_in_set_constraints = true))
#     println("number of constraints: $(JuMP.num_constraints(intersecting_polytopes_jump[i]; count_variable_in_set_constraints = true))")
#     println()   
#     println()
#     readline()
# end
# TODO: 
# 1) implement deeper intersection, by moving polytopes towards the average of convex hull vertices? Maybe use https://github.com/JuliaPolyhedra/Polyhedra.jl/blob/8c131a6cc883c541922a8b8efe835932e8b2593f/src/center.jl#L7
# 2) also, implement distance function calculator, so we know the optimal solution

# - USE JUMP MODELS OF POLYTOPES INTO FRANKWOLFE.JL LMOS
# - TRY IF EVERYTHING WORKS WITH BPCG