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
function generate_polytope(config::Config, n_points::Int, bounds::Vector{Tuple{Float64, Float64}})
    # Check that the number of columns specified in `config.n` matches the length of `bounds`
    if config.n != length(bounds)
        error("Mismatch: The number of columns in `config` (config.n = $(config.n)) must match the number of elements in `bounds` (length = $(length(bounds))).")
    end

    # Initialize a matrix of zeros with the given number of points (rows) and columns
    vertices = zeros(Float64, n_points, config.n)

    # Fill each column based on the specified bounds
    for j in 1:config.n
        lower_bound, upper_bound = bounds[j]
        # Generate random numbers within the given bounds
        vertices[:, j] = lower_bound .+ (upper_bound - lower_bound) .* rand(Float64, n_points)
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

    return _, v2_closest, min_dist
end

# Function to intersect two polytopes, moving the second onto the closest vertex of the first
function intersect_polytopes(config::Config, V1::Matrix{T}, V2::Matrix{T}, polytope_to_move::Polyhedron{T}) where T
    
    # Find the closest pair of points between the two sets of points
    v1_closest, v2_closest, _ = closest_pair(config, V1, V2)

    # Create a translation vector (offset) from p_closest - q_closest
    direction = v1_closest - v2_closest

    # Create "fake" offset polytope
    O = polyhedron(vrep([direction]), CDDLib.Library())

    # Return translated polytope
    return polytope_to_move + O
end
# (Multiple Dispatch) Function to intersect two polytopes, moving the second unto a given vertex of the first
function intersect_polytopes(config::Config, v::Vector, V2::Matrix{T}, polytope_to_move::Polyhedron{T}) where T
    
    # Find the closest pair of points between the two sets of points
    _, v2_closest, _ = closest_pair(config, v, V2)

    # Create a translation vector (offset) from p_closest - q_closest
    direction = v - v2_closest

    # Create "fake" offset polytope
    O = polyhedron(vrep([direction]), CDDLib.Library())

    # Return translated polytope
    return polytope_to_move + O
end

# Function to plot polytopes with random colors
function plot_polytopes(polytopes::Vector{Polyhedra.Polyhedron})
    for poly in polytopes
        # Generate a random RGB color
        random_color = RGB(rand(), rand(), rand())

        # Plot the polytope with the random color
        plot!(poly, color=random_color, alpha=0.6)
    end
end

# Function to generate a polytope with a given JuMP model
function polyhedra_to_jump(config::Config, polytope::Polyhedra.Polyhedron{T}) where T
    model = Model()
    @variable(model, x[1:config.n])
    @constraint(model, x in polytope)
    
    return model
end

"""
TODO: 
1) CHANGE INTERSECTION: THE THING YOU'RE DOING RIGHT NOW DOESN'T WORK FOR K >= 3 --> check"!
2) use planar_hull and or remove redundancy!!!
2) USE JUMP MODELS OF POLYTOPES INTO FRANKWOLFE.JL LMOS
3) MAYBE GENERATE POLYTOPE SHOULD BE IN /SRC?
4) MAYBE DECIDE NUMBER OF VERTICES OF EACH POLYTOPE IN CONFIG? OR GENERATE IT RANDOMLY AS SOMETHING >= N+1 INSIDE POLYTOPES.JL
4) TRY IF EVERYTHING WORKS WITH BPCG
"""