# `polytopes.jl`
# Function to generate non-overlapping bounds for multiple dimensions within [1e-03, 1e+03]
# `margin`: ensures min distance between current polytope's UB and next polytope's LB, to create non-overlapping polytopes
function generate_nonintersecting_bounds(config::Config; margin::Float64=10.0, stepsize::Float64=200.0, start_point::Float64=-100.0)
    
    bounds_list = Vector{Vector{Tuple{Float64, Float64}}}(undef, config.k)

    # Generate the first bounds randomly as pairs (LB, UB), within given range
    lower_bound = rand(start_point : start_point + stepsize - margin)
    upper_bound = rand(lower_bound + margin : lower_bound + stepsize + margin)
    bounds_list[1] = [(lower_bound, upper_bound) for _ in 1:config.n]

    # Generate subsequent bounds while ensuring no overlap
    for i in 2:config.k
        bounds = Vector{Tuple{Float64, Float64}}(undef, config.n)

        # Initialize bounds for current polytope
        for d in 1:config.n
            # Check bounds from previous polytope
            _, prev_max = bounds_list[i - 1][d]
            # Ensure the new lower bound is at least 'margin' away from the previous upper bound
            lower_bound = rand(prev_max + margin : prev_max + stepsize + margin)
            # Generate the upper bound ensuring it stays within the range
            upper_bound = rand(lower_bound + margin : lower_bound + stepsize + margin)
            # Assign new bounds for current dimension
            bounds[d] = (lower_bound, upper_bound)
        end
        # Add bounds for current polytope
        bounds_list[i] = bounds
    end

    return bounds_list
end

# TODO: THIS CAUSES NUMERICAL ISSUES WITH THE Polyhedra.jl LIBRARY, DON'T USE FOR NOW (relevant functions commented)
# Create Polyhedra.Polyhedron from list of vertices
function polytope(vertices::Matrix{T}) where T
    
    # return polyhedron(vrep(vertices), CDDLib.Library())
    return polyhedron(vrep(vertices), CDDLib.Library(:exact))
    # return polyhedron(vrep(vertices))
end
# (Multiple dispatch)
function polytope(vertices::Vector{Vector{T}}) where T
    
    # return polyhedron(vrep(vertices), CDDLib.Library())
    return polyhedron(vrep(vertices), CDDLib.Library(:exact))
    # return polyhedron(vrep(vertices))
end
# TODO: THIS CAUSES NUMERICAL ISSUES WITH THE Polyhedra.jl LIBRARY, DON'T USE FOR NOW (relevant functions commented)
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
    
    return vertices
end

# Generate one polytope: generate one random set of vertices, within given bounds for each dimension
function generate_polytope(config::Config, idx::Int, bounds::Vector{Tuple{T, T}}) where T
    
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

    # TODO: THIS CAUSES NUMERICAL ISSUES WITH THE Polyhedra.jl LIBRARY, DON'T USE FOR NOW (relevant functions commented)
    # Generate polytope as convex hull of given points, and remove redundant (i.e. non-vertex) points
    #vertices = nonredundant_polytope(vertices)
    
    return vertices
end

# Generate k nonintersecting polytope: generate k random sets of vertices, within given bounds for each dimension
function generate_nonintersecting_polytopes(config::Config, bounds::Vector{Vector{Tuple{T, T}}}) where T
    
    vertices = Vector{Matrix{Float64}}()

    # Generate k polytopes within the specified bounds
    for i in 1:config.k
        # Generate vertices and corresponding polytope
        verts = generate_polytope(config, i, bounds[i])

        # Append to `vertices` and `polytopes`
        push!(vertices, verts)
    end

    nonintersecting_polytopes_lmos = create_lmos(config, vertices)
    _, _, primal, fw_gap = compute_distance(config, nonintersecting_polytopes_lmos)

    # Check that polytope intersection is empty
    if approxequal(primal, 0.0)
        error("Invalid polytopes: they should not intersect, but the total distance among them is $primal.")
    end

    return vertices, primal, fw_gap
end

# Generate nonintersecting polytopes, then intersect them in (at least) one point, then compute total distance between them
function generate_polytopes(config::Config)

    # `n_points` contains n. of vertices used to generate each polytope
    println("Generating $(config.k) intersecting polytopes of dimension $(config.n) (n. vertices for each: $(config.n_points))")
    
    # Generate random, non-intersecting bounds
    bounds_list = generate_nonintersecting_bounds(config)
    
    # Generate non intersecting polytopes
    vertices, primal, fw_gap = generate_nonintersecting_polytopes(config, bounds_list)
        
    # Generate intersecting polytopes
    shifted_vertices = intersect_polytopes(config, vertices)
    
    # Check that polytope intersection is not empty
    lmos_shifted = create_lmos(config, shifted_vertices)
    _, _, primal_shifted, _ = compute_distance(config, lmos_shifted)
    
    if !approxequal(primal_shifted, 0.0)
        error("Invalid polytopes: they should intersect, but the total distance among them is $primal_shifted")
    end

    return vertices, shifted_vertices, primal, fw_gap
end

# Function to find the closest points between two sets of points
function closest_pair(config::Config, V1::Matrix{T}, V2::Matrix{T}) where T

    # Check that both matrices have the correct number of columns
    if size(V1, 2) != config.n || size(V2, 2) != config.n
        error("Both point sets must have the dimension specified in `config.n`.")
    end
    
    lmos = create_lmos(config, [V1, V2])
    _, last_lmo_solution, _, _ = compute_distance(config, lmos)
    # `last_lmo_solution` is a FrankWolfe.BlockVector
    v1_closest, v2_closest = last_lmo_solution.blocks[1], last_lmo_solution.blocks[2]

    return v1_closest, v2_closest
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

    return v2_closest
end

# Function to intersect two polytopes, moving the second onto the closest vertex of the first
function intersect_polytopes(config::Config, V1::Matrix{T}, V2::Matrix{T}) where T
    
    # Find the closest pair of points between the two sets of points
    v1_closest, v2_closest = closest_pair(config, V1, V2)

    # Create offset vector `v1_closest - v2_closest`
    direction = v1_closest - v2_closest

    # Shift each vertex in V2 along `direction`
    shifted_vertices_curr = V2 .+ direction'

    # Return translated polytope and other info (`points()` returns a matrix with the vertices)
    return shifted_vertices_curr, v1_closest, v2_closest
end
# (Multiple Dispatch) Function to intersect two polytopes, moving the second onto a given vertex of the first
function intersect_polytopes(config::Config, v::Vector{T}, V2::Matrix{T}) where T
    
    # Find the closest pair of points between the two sets of points
    v2_closest = closest_pair(config, v, V2)

    # Create offset vector `v - v2_closest`
    direction = v - v2_closest

    # Shift each vertex in V2 along `direction`
    shifted_vertices_curr = V2 .+ direction'

    # Return translated polytope and other info (`points()` returns a matrix with the vertices)
    return shifted_vertices_curr, v2_closest
end
# (Multiple Dispatch) Function to intersect k polytopes so that the intersection contains (at least) one point
function intersect_polytopes(
    config::Config,
    vertices::Vector{Matrix{T}},
    ) where T
    
    # Initialize empty vectors for vertices, Polyhedra polytopes, and Polyhedra intersecting polytopes, JuMP intersecting polytopes
    shifted_vertices = Vector{Matrix{Float64}}()

    # Update P₁ data
    push!(shifted_vertices, vertices[1])
    
    # Move P₂ towards P₁ so that they intersect (at least) in v₁, then update P₂ data
    shifted_vertices_curr, v₁, _ = intersect_polytopes(config, vertices[1], vertices[2])
    push!(shifted_vertices, shifted_vertices_curr)

    # Move each subsequent polytope Pᵢ, for k > 3, to intersect with P₁ in at least v₁, then update Pᵢ data
    for i in 3:config.k
        # Move iᵗʰ polytope so that it intersects with the others in (at least) v₁
        shifted_vertices_curr, _ = intersect_polytopes(config, v₁, vertices[i])

        # Update data
        push!(shifted_vertices, shifted_vertices_curr)
    end
    # check_intersection(intersecting_polytopes_polyhedra)

    return shifted_vertices
    
end

# Function to generate a JuMP.Model object from a Polyhedra.Polyhedron object
function polyhedra_to_jump(config::Config, polytope::Polyhedron{T}) where T
    
    model = Model(SCIP.Optimizer)
    set_optimizer_attribute(model, "display/verblevel", 0)
    @variable(model, x[1:config.n])
    @constraint(model, x in polytope)
    # Ensure the model is optimized (needed for reusability of the model in FrankWolfe algorithms)
    optimize!(model)  

    return model
end
# (Multiple dispatch) Function to generate a Vector{JuMP.Model} from a Vector{Polyhedra.Polyhedron}
function polyhedra_to_jump(config::Config, polytopes::Vector{Polyhedron{T}}) where T
    
    polytopes_jump = Vector{Model()}()
    for polytope in polytopes
        polytope_jump = polyhedra_to_jump(config, polytope)
        push!(polytopes_jump, polytope_jump)
    end
    return polytopes_jump
end

# TODO: THIS CAUSES NUMERICAL ISSUES WITH THE Polyhedra.jl LIBRARY, DON'T USE FOR NOW (relevant functions commented)
# Check that intersection of `intersecting_polytopes` is not empty
function check_intersection(intersecting_polytopes::Vector{Polyhedron})
    
    # Get the number of points in the intersection of all polytopes
    intersection_size = npoints(intersect(intersecting_polytopes...))
    
    # Print intersection size for debugging
    println("Intersection of all polytopes is a polytope with $intersection_size vertices")
    
    # Assert that the intersection is not empty
    @assert intersection_size ≥ 1 "There must be at least one point in the intersection of all polytopes"
end

# Compute distance between k polytopes, by running the FW algorithm
function compute_distance(config::Config, lmo_list::Vector{FrankWolfe.LinearMinimizationOracle})

    # Redefine `config.max_iterations`
    config_opt = Config("examples/config.yml"; max_iterations=config.max_iterations_opt)

    # Create FrankWolfe.ProductLMO from list of LMOs
    prod_lmo = create_product_lmo(lmo_list)
    
    # Run Block-coordinate BPCG with CyclicUpdate
    last_iterate, last_lmo_vertex, primal, fw_gap, _ = run_FW(config_opt, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
    
    return last_iterate, last_lmo_vertex, primal, fw_gap
end
# (Multiple dispatch) Compute distance between k polytopes, by running the FW algorithm
function compute_distance(config::Config, lmo_list::Vector{FrankWolfe.ConvexHullOracle})

    # Redefine `config.max_iterations`
    config_opt = Config("examples/config.yml"; max_iterations=config.max_iterations_opt)
    
    # Create FrankWolfe.ProductLMO from list of LMOs
    prod_lmo = create_product_lmo(lmo_list)
    
    # Run Block-coordinate BPCG with CyclicUpdate
    last_iterate, last_lmo_vertex, primal, fw_gap, _ = run_FW(config_opt, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
    
    return last_iterate, last_lmo_vertex, primal, fw_gap
end
# (Multiple dispatch) Compute distance between k polytopes, by running the FW algorithm
function compute_distance(config::Config, lmo_list::Vector{FrankWolfe.MathOptLMO})

    # Redefine `config.max_iterations`
    config_opt = Config("examples/config.yml"; max_iterations=config.max_iterations_opt)

    # Create FrankWolfe.ProductLMO from list of LMOs
    prod_lmo = create_product_lmo(lmo_list)
    
    # Run Block-coordinate BPCG with CyclicUpdate
    last_iterate, last_lmo_vertex, primal, fw_gap, _ = run_FW(config_opt, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
    
    return last_iterate, last_lmo_vertex, primal, fw_gap
end

# Save data to given .jld2 file
function save_polytopes(
    filename::String,
    vertices::Vector{Matrix{T}},
    shifted_vertices::Vector{Matrix{T}},
    primal::T,
    fw_gap::T
    ) where T

    save(filename, Dict("vertices" => vertices, "shifted_vertices" => shifted_vertices, "primal" => primal, "fw_gap" => fw_gap))
    println("Saving data to $filename")
end
# (Multiple dispatch) Automatically generated .jld2 filename
function save_polytopes(
    config::Config,
    vertices::Vector{Matrix{T}},
    shifted_vertices::Vector{Matrix{T}},
    primal::T,
    fw_gap::T
    ) where T
    
    filename = generate_filename(config, vertices)
    save(filename, Dict("vertices" => vertices, "shifted_vertices" => shifted_vertices, "primal" => primal, "fw_gap" => fw_gap))
    println("Saving data to $filename")
    
    return filename
end

# Load data from .jld2 file
function load_polytopes(filename::String)
    
    f = load(filename)
    vertices = f["vertices"]
    shifted_vertices = f["shifted_vertices"]
    primal = f["primal"]
    fw_gap = f["fw_gap"]
    
    return vertices, shifted_vertices, primal, fw_gap
end