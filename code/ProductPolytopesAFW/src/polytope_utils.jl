# `polytope_utils.jl`

# TODO: THIS CAUSES NUMERICAL ISSUES WITH THE Polyhedra.jl LIBRARY, DON'T USE FOR NOW (relevant functions commented)
# Create Polyhedra.Polyhedron from list of vertices
function polytope(vertices::Matrix{T}) where T
    
    # return polyhedron(vrep(vertices), CDDLib.Library())
    return polyhedron(vrep(vertices), CDDLib.Library(:exact))
    # return polyhedron(vrep(vertices))
end

# (Multiple Dispatch) Create Polyhedra.Polyhedron from list of vertices
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

# Function to compute the vertex mean (average of sampled vertices) of a given list of vertices
function vertex_mean(vertices::Matrix{T}) where T
    # Check if the matrix is empty
    if size(vertices, 1) == 0
        error("The matrix of vertices is empty.")
    end
    
    _vertex_mean = Vector{T}()

    # Sum all the vertices
    sum_vector = vec(sum(vertices, dims=1))
    
    # Compute the vertex mean by dividing by the number of vertices
    _vertex_mean = sum_vector ./ size(vertices, 1)
    
    return _vertex_mean
end

# Sample a random vertex (row) from a vertex matrix (set of vertices)
function random_vertex(vertices::Matrix{T}, rng::AbstractRNG) where T
    nverts = size(vertices, 1)
    if nverts == 0
        error("Cannot sample a random vertex: the vertex matrix is empty.")
    end
    idx = rand(rng, 1:nverts)
    return vec(vertices[idx, :])
end

# Sample a random convex combination of the given vertices (uniform over the simplex weights)
function random_convex_combination(vertices::Matrix{T}, rng::AbstractRNG) where T
    nverts = size(vertices, 1)
    if nverts == 0
        error("Cannot sample a random convex combination: the vertex matrix is empty.")
    end

    weights = rand(rng, nverts)
    s = sum(weights)
    if s == 0
        weights .= 1 / nverts # extremely unlikely --> fall back to uniform weights
    else
        weights ./= s # normalize
    end

    combo = weights' * vertices # 1xn
    return vec(combo)
end

# Compute the mean of per-polytope vertex means in `vertices`
function global_mean_of_centers(config::Config, vertices::Vector{Matrix{T}}) where T
    if isempty(vertices)
        error("Cannot compute global mean of centers: no polytopes provided.")
    end
    centers = [vertex_mean(V) for V in vertices]
    A = reduce(hcat, centers) # stack k centers as columns into an nxk matrix
    return vec(mean(A, dims=2)) # n-element vector
end

# (Multiple Dispatch) Find closest points between two polytopes by running FW
function closest_pair(config::Config, V1::Matrix{T}, V2::Matrix{T}) where T

    # Check that both matrices have the correct number of columns
    if size(V1, 2) != config.n || size(V2, 2) != config.n
        error("Both point sets must have the dimension specified in `config.n`.")
    end
    
    lmos = create_lmos(config, [V1, V2])
    _, last_lmo_solution, _, _ = compute_distance(config, lmos)
    # `last_lmo_solution` is a FrankWolfe.BlockVector
    # Materialize as dense vectors: with matrix-backed ConvexHullLMO, blocks may be views.
    v1_closest = Vector{T}(last_lmo_solution.blocks[1])
    v2_closest = Vector{T}(last_lmo_solution.blocks[2])

    return v1_closest, v2_closest
end

# (Multiple Dispatch) Compute a reference point for a single polytope given an anchor
function compute_reference_point(config::Config, anchor::AbstractVector{T}, vertices::Matrix{T}) where T
    if config.intersection_reference_point == "center"
        return vertex_mean(vertices) # ignores `anchor`
    elseif config.intersection_reference_point == "vertex"
        # Closest vertex of this polytope to the anchor
        return closest_pair(config, anchor, vertices)
    else
        error("Unknown intersection reference point: $(config.intersection_reference_point)")
    end
end
# (Multiple Dispatch) Compute reference points for all polytopes given an anchor
function compute_reference_points(config::Config, vertices::Vector{Matrix{T}}, anchor::AbstractVector{T}) where T
    nP = length(vertices)
    refs = Vector{Vector{T}}(undef, nP) # preallocate one reference point per polytope
    for i in 1:nP
        refs[i] = compute_reference_point(config, anchor, vertices[i])
    end
    return refs
end

# Shift a polytope so that `reference` moves towards `anchor` by a given stepsize
function shift_polytope(
    vertices::Matrix{T},
    reference::AbstractVector{T},
    anchor::AbstractVector{T};
    stepsize::T = one(T),
) where T
    if size(vertices, 2) != length(reference) || length(reference) != length(anchor)
        error("Dimension mismatch: vertices have $(size(vertices, 2)) columns, but reference/anchor have lengths $(length(reference)) and $(length(anchor)).")
    end
    # Translation direction; with stepsize = 1, `anchor` lies in the shifted polytope
    direction = stepsize * (anchor .- reference)
    return vertices .+ direction'
end

# (Multiple Dispatch) Here the first input is a vector (single point)
function closest_pair(config::Config, v::AbstractVector{T}, V2::Matrix{T}) where T
    
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
    v2_closest = Vector{T}(V2[closest_index, :])

    return v2_closest
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

# (Multiple Dispatch) Function to generate a Vector{JuMP.Model} from a Vector{Polyhedra.Polyhedron}
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

# (Multiple Dispatch) Compute distance between k polytopes, by running the FW algorithm
function compute_distance(config::Config, lmo_list::Vector{FrankWolfe.LinearMinimizationOracle})

    # Redefine `config.max_iterations`
    config_opt = modify_config(config, max_iterations=config.max_iterations_opt)

    # Create FrankWolfe.ProductLMO from list of LMOs
    prod_lmo = create_product_lmo(lmo_list)

    # Run Block-coordinate BPCG with CyclicUpdate
    # last_iterate, last_lmo_vertex, primal, fw_gap, _ = run_BlockCoordinateFW(config_opt, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
    last_iterate, last_lmo_vertex, primal, fw_gap, _ = run_BlockCoordinateFW(config_opt, FrankWolfe.CyclicUpdate(), AwayStep(), prod_lmo)
    
    return last_iterate, last_lmo_vertex, primal, fw_gap
end
# (Multiple Dispatch) Compute distance between k polytopes, by running the FW algorithm
function compute_distance(config::Config, lmo_list::Vector{FrankWolfe.ConvexHullLMO})

    # Redefine `config.max_iterations`
    config_opt = modify_config(config, max_iterations=config.max_iterations_opt)

    # Create FrankWolfe.ProductLMO from list of LMOs
    prod_lmo = create_product_lmo(lmo_list)

    # Run Block-coordinate BPCG with CyclicUpdate
    # last_iterate, last_lmo_vertex, primal, fw_gap, _ = run_BlockCoordinateFW(config_opt, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
    last_iterate, last_lmo_vertex, primal, fw_gap, _ = run_BlockCoordinateFW(config_opt, FrankWolfe.CyclicUpdate(), AwayStep(), prod_lmo)


    return last_iterate, last_lmo_vertex, primal, fw_gap
end
# (Multiple Dispatch) Compute distance between k polytopes, by running the FW algorithm
function compute_distance(config::Config, lmo_list::Vector{FrankWolfe.MathOptLMO})

    # Redefine `config.max_iterations`
    config_opt = modify_config(config, max_iterations=config.max_iterations_opt)

    # Create FrankWolfe.ProductLMO from list of LMOs
    prod_lmo = create_product_lmo(lmo_list)

    # Run Block-coordinate BPCG with CyclicUpdate
    # last_iterate, last_lmo_vertex, primal, fw_gap, _ = run_BlockCoordinateFW(config_opt, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
    last_iterate, last_lmo_vertex, primal, fw_gap, _ = run_BlockCoordinateFW(config_opt, FrankWolfe.CyclicUpdate(), AwayStep(), prod_lmo)

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
    
    filename = generate_filename(config)
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
