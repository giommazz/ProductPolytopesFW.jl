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

# Create Polyhedra.Polyhedron from list of vertices
function polytope(vertices::Matrix{T}) where T
    return polyhedron(vrep(vertices), CDDLib.Library())
end
# (Multiple dispatch)
function polytope(vertices::Vector{Vector{T}}) where T
    return polyhedron(vrep(vertices), CDDLib.Library())
end

# Create non-redundant Polyhedra.Polyhedron from list of vertices
function nonredundant_polytope(vertices::Matrix{T}; redundancy_flag=true::Bool) where T
    # Use vertices to create object of type Polyhedra.Polyhedron 
    poly = polytope(vertices)
    if redundancy_flag
        # Find vertices of `poly`
        removevredundancy!(poly)
    end
    # Type conversion to obtain Matrix{T} from `Polyhedra.points(poly)`
    vertices = stack(collect(points(poly)), dims=1)
    
    return vertices, poly
end

# Function to generate a random vertex matrix with bounds for each column
function generate_polytopes(config::Config, idx::Int, bounds::Vector{Tuple{T, T}}) where T
    # Check that the number of columns specified in `config.n` matches the length of `bounds`
    if config.n != length(bounds)
        error("Mismatch: The number of columns in `config` (config.n = $(config.n)) must match the number of elements in `bounds` (length = $(length(bounds))).")
    end

    # Initialize a matrix of zeros with the given number of points (rows) and columns
    vertices = zeros(Float64, config.n_points[idx], config.n)

    # Fill each column based on the specified bounds
    for j in 1:config.n
        lower_bound, upper_bound = bounds[j]
        # Generate vector containing `config.n_points[idx]` random Float64 within the given bounds
        vertices[:, j] = lower_bound .+ (upper_bound - lower_bound) .* rand(Float64, config.n_points[idx])
    end
    # Generate polytope as convex hull of given points, and remove redundant (i.e. non-vertex) points
    vertices, poly = nonredundant_polytope(vertices)
    
    return vertices, poly
end
# Function to generate k random vertex matrices with bounds for each column
function generate_polytopes(config::Config, bounds::Vector{Vector{Tuple{T, T}}}, vertices::Vector{Matrix{T}}, polytopes::Vector{Polyhedron}) where T
    # Generate k polytopes within the specified bounds
    for i in 1:config.k
        # Generate vertices and corresponding polytope
        verts, poly = generate_polytopes(config, i, bounds[i])
        
        # Append to `vertices` and `polytopes`
        push!(vertices, verts)
        push!(polytopes, poly)
    end
end

# Function to find the closest points between two sets of points
function closest_pair(config::Config, V1::Matrix{T}, V2::Matrix{T}) where T

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
function closest_pair(config::Config, v::Vector{T}, V2::Matrix{T}) where T
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
    P = polytope([direction])
    # Shift `polytope_to_move` along the offset, to obtain a shifted polytope
    shifted_polytope_curr = polytope_to_move + P
    # Apply the translation to each vertex in V2 to get shifted_vertices
    shifted_vertices_curr = V2 .+ direction'

    # Return translated polytope and other info (`points()` returns a matrix with the vertices)
    return shifted_polytope_curr, shifted_vertices_curr, v1_closest, v2_closest, distance
end
# (Multiple Dispatch) Function to intersect two polytopes, moving the second unto a given vertex of the first
function intersect_polytopes(config::Config, v::Vector{T}, V2::Matrix{T}, polytope_to_move::Polyhedron{T}) where T
    
    # Find the closest pair of points between the two sets of points
    v2_closest, distance = closest_pair(config, v, V2)

    # Create a translation vector (offset) from p_closest - q_closest
    direction = v - v2_closest

    # Create "fake" offset polytope
    P = polytope([direction])
    # Shift `polytope_to_move` along the offset
    shifted_polytope_curr = polytope_to_move + P
    # Apply the translation to each vertex in V2 to get shifted_vertices
    shifted_vertices_curr = V2 .+ direction'

    # Return translated polytope and other info (`points()` returns a matrix with the vertices)
    return shifted_polytope_curr, shifted_vertices_curr, v2_closest, distance
end

# (Multiple Dispatch) Function to intersect k polytopes
function intersect_polytopes(
    config::Config,
    vertices::Vector{Matrix{T}},
    shifted_vertices::Vector{Matrix{T}},
    polytopes::Vector{Polyhedron},
    intersecting_polytopes_polyhedra::Vector{Polyhedron},
    intersecting_polytopes_jump::Vector{Model}
    ) where T
    
    # Update P₁ data
    push!(intersecting_polytopes_polyhedra, polytopes[1])
    push!(intersecting_polytopes_jump, polyhedra_to_jump(config, polytopes[1]))
    push!(shifted_vertices, vertices[1])
    
    # Move P₂ towards P₁ so that they intersect (at least) in v₁, then update P₂ data
    shifted_polytope_curr, shifted_vertices_curr, v₁, _, distance = intersect_polytopes(config, vertices[1], vertices[2], polytopes[2])
    push!(intersecting_polytopes_polyhedra, shifted_polytope_curr)
    push!(intersecting_polytopes_jump, polyhedra_to_jump(config, shifted_polytope_curr))
    push!(shifted_vertices, shifted_vertices_curr)

    # Move each subsequent polytope Pᵢ, for k > 3, to intersect with P₁ in at least v₁, then update Pᵢ data
    for i in 3:config.k
        # Move iᵗʰ polytope so that it intersects with the others in (at least) v₁
        shifted_polytope_curr, shifted_vertices_curr, _, distance = intersect_polytopes(config, v₁, vertices[i], polytopes[i])

        # Update data
        push!(intersecting_polytopes_polyhedra, shifted_polytope_curr)
        push!(intersecting_polytopes_jump, polyhedra_to_jump(config, shifted_polytope_curr))
        push!(shifted_vertices, shifted_vertices_curr)
    end

    #check_intersection(intersecting_polytopes_polyhedra)
end

# Function to generate a polytope with a given JuMP model
function polyhedra_to_jump(config::Config, polytope::Polyhedron{T}) where T

    model = Model(GLPK.Optimizer)
    @variable(model, x[1:config.n])
    @constraint(model, x in polytope)
    # Ensure the model is optimized (neede for reusability of the model)
    optimize!(model)  
    return model
end

# Check that intersection of `intersecting_polytopes` is not empty
# Check that intersection of `intersecting_polytopes` is not empty
function check_intersection(intersecting_polytopes::Vector{Polyhedron})
    
    # Get the number of points in the intersection of all polytopes
    intersection_size = npoints(intersect(intersecting_polytopes...))
    
    # Print intersection size for debugging
    println("Intersection of all polytopes is a polytope with $intersection_size vertices")
    
    # Assert that the intersection is not empty
    @assert intersection_size ≥ 1 "There must be at least one point in the intersection of all polytopes"
end

# Main function to generate and move polytopes
# `n_points` contains n. of vertices used to generate each polytope
function generate_intersecting_polytopes(config::Config)

    println("Generating $(config.k) intersecting polytopes of dimension $(config.n) (n. vertices for each: $(config.n_points))")
    
    # Generate random, non-intersecting bounds
    bounds_list = generate_non_intersecting_bounds(config)

    # Initialize empty vectors for vertices, Polyhedra polytopes, and Polyhedra intersecting polytopes, JuMP intersecting polytopes
    vertices = Vector{Matrix{Float64}}()
    shifted_vertices = Vector{Matrix{Float64}}()
    polytopes = Vector{Polyhedron}()
    intersecting_polytopes_polyhedra = Vector{Polyhedron}()
    intersecting_polytopes_jump = Vector{Model}()

    # Generate non intersecting polytopes
    generate_polytopes(config, bounds_list, vertices, polytopes)

    intersect_polytopes(config, vertices, shifted_vertices, polytopes, intersecting_polytopes_polyhedra, intersecting_polytopes_jump)

    return vertices, shifted_vertices, polytopes, intersecting_polytopes_polyhedra, intersecting_polytopes_jump
end

# Save data to given .jld2 file
function save_intersecting_polytopes(
    filename::String,
    vertices::Vector{Matrix{T}},
    shifted_vertices::Vector{Matrix{T}}
    ) where T

    # TODO: SALVA VERTICES E VERTICES_INTERSECTING: GABA
    save(filename, Dict("vertices" => vertices, "shifted_vertices" => shifted_vertices))
    println("Saving data to $filename")
end
# (Multiple dispatch) Automatically generated .jld2 filename
function save_intersecting_polytopes(
    config::Config,
    vertices::Vector{Matrix{T}},
    shifted_vertices::Vector{Matrix{T}}
    ) where T
    
    filename = generate_filename(config, vertices)
    save(filename, Dict("vertices" => vertices, "shifted_vertices" => shifted_vertices))
    println("Saving data to $filename")
    return filename
end

# Load data from .jld2 file
function load_intersecting_polytopes(filename::String)
    f = load(filename)
    vertices = f["vertices"]
    shifted_vertices = f["shifted_vertices"]
    return vertices, shifted_vertices
end

# for i in 1:config.k
#     println("°°°°°°°°°°°°°°°°---------------------------------> [$i]")
    
#     println("\tVREP")
#     for constraint in list_of_constraint_types(intersecting_polytopes_jump[i])
#         println("\t\tConstraint: $constraint")
#     end
#     println()
#     println("\t\tall constraints: ", all_constraints(intersecting_polytopes_jump[i]; include_variable_in_set_constraints = true))
#     println("number of constraints: $(num_constraints(intersecting_polytopes_jump[i]; count_variable_in_set_constraints = true))")
#     println()   
#     println()


#     PP = hrep(intersecting_polytopes_polyhedra[i])
#     println("\tHREP")
#     PPJ = polyhedra_to_jump(config, intersecting_polytopes_polyhedra[i])
#     for constraint in list_of_constraint_types(PPJ)
#         println("\t\tConstraint: $constraint")
#     end
#     println()
#     println("\t\tall constraints: ", all_constraints(intersecting_polytopes_jump[i]; include_variable_in_set_constraints = true))
#     println("number of constraints: $(num_constraints(intersecting_polytopes_jump[i]; count_variable_in_set_constraints = true))")
#     println()   
#     println()
#     readline()
# end
# TODO: 
# 1) implement deeper intersection, by moving polytopes towards the average of convex hull vertices? Maybe use https://github.com/JuliaPolyhedra/jl/blob/8c131a6cc883c541922a8b8efe835932e8b2593f/src/center.jl#L7
# 2) also, implement distance function calculator, so we know the optimal solution