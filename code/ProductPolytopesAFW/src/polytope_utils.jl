# `polytope_utils.jl`

"""
    polytope(vertices)

Create a `Polyhedra.Polyhedron` from a list/matrix of vertices (V-representation).

Note: this path has caused numerical issues with `Polyhedra.jl`/`CDDLib` in the past; use with care.
"""
function polytope(vertices::Matrix{T}) where T
    
    # return polyhedron(vrep(vertices), CDDLib.Library())
    return polyhedron(vrep(vertices), CDDLib.Library(:exact))
    # return polyhedron(vrep(vertices))
end

"""
    polytope(vertices)

Multiple dispatch: build a `Polyhedra.Polyhedron` from vertices provided as a `Vector{Vector}`.
"""
function polytope(vertices::Vector{Vector{T}}) where T
    
    # return polyhedron(vrep(vertices), CDDLib.Library())
    return polyhedron(vrep(vertices), CDDLib.Library(:exact))
    # return polyhedron(vrep(vertices))
end

"""
    nonredundant_polytope(vertices; redundancy_flag=true) -> Matrix

Create a `Polyhedra.Polyhedron` from `vertices` and (optionally) remove redundant points, returning the
remaining vertices as a matrix.

Note: this path has caused numerical issues with `Polyhedra.jl`/`CDDLib` in the past; use with care.
"""
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

"""
    vertex_mean(vertices) -> Vector

Compute the mean of the rows of `vertices` (average sampled vertex).
"""
function vertex_mean(vertices::Matrix{T}) where T
    # Check if the matrix is empty
    if size(vertices, 1) == 0
        error("The matrix of vertices is empty.")
    end
    
    _vertex_mean = Vector{T}()

    # Sum all vertices
    sum_vector = vec(sum(vertices, dims=1))
    
    # Compute the vertex mean by dividing by the number of vertices
    _vertex_mean = sum_vector ./ size(vertices, 1)
    
    return _vertex_mean
end

"""
    random_vertex(vertices, rng) -> Vector

Sample one random row from `vertices` using `rng`.
"""
function random_vertex(vertices::Matrix{T}, rng::AbstractRNG) where T
    nverts = size(vertices, 1)
    if nverts == 0
        error("Cannot sample a random vertex: the vertex matrix is empty.")
    end
    idx = rand(rng, 1:nverts)
    return vec(vertices[idx, :])
end

"""
    random_convex_combination(vertices, rng) -> Vector

Sample a random convex combination of the rows of `vertices` (weights uniform on the simplex).
"""
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

"""
    global_mean_of_centers(config, vertices) -> Vector

Compute the mean of per-polytope vertex means across a list of vertex matrices.
"""
function global_mean_of_centers(config::Config, vertices::Vector{Matrix{T}}) where T
    if isempty(vertices)
        error("Cannot compute global mean of centers: no polytopes provided.")
    end
    centers = [vertex_mean(V) for V in vertices]
    A = reduce(hcat, centers) # stack k centers as columns into an nxk matrix
    return vec(mean(A, dims=2)) # n-element vector
end

"""
    closest_pair(config, V1, V2) -> (v1_closest, v2_closest)

Multiple dispatch: find a closest pair of points between two polytopes (given by their vertices),
by solving the distance problem with FW.
"""
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

"""
    compute_reference_point(config, anchor, vertices) -> Vector

Multiple dispatch: compute the reference point for one polytope given the chosen `anchor`.
"""
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
"""
    compute_reference_points(config, vertices, anchor) -> Vector{Vector}

Multiple dispatch: compute one reference point per polytope given the chosen `anchor`.
"""
function compute_reference_points(config::Config, vertices::Vector{Matrix{T}}, anchor::AbstractVector{T}) where T
    nP = length(vertices)
    refs = Vector{Vector{T}}(undef, nP) # preallocate one reference point per polytope
    for i in 1:nP
        refs[i] = compute_reference_point(config, anchor, vertices[i])
    end
    return refs
end

"""
    shift_polytope(vertices, reference, anchor; stepsize=1) -> Matrix

Translate `vertices` so that `reference` moves towards `anchor` by the given `stepsize`.
"""
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

"""
    closest_pair(config, v, V2) -> Vector

Multiple dispatch: return the point in `V2` closest to `v` (brute-force scan).
"""
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

"""
    polyhedra_to_jump(config, polytope) -> JuMP.Model

Convert a `Polyhedra.Polyhedron` into a `JuMP.Model` with constraint `x in polytope`.
"""
function polyhedra_to_jump(config::Config, polytope::Polyhedron{T}) where T
    
    model = Model(SCIP.Optimizer)
    set_optimizer_attribute(model, "display/verblevel", 0)
    @variable(model, x[1:config.n])
    @constraint(model, x in polytope)
    # Ensure the model is optimized (needed for reusability of the model in FrankWolfe algorithms)
    optimize!(model)  

    return model
end

"""
    polyhedra_to_jump(config, polytopes) -> Vector{JuMP.Model}

Multiple dispatch: convert a list of `Polyhedra.Polyhedron` objects into `JuMP.Model`s.
"""
function polyhedra_to_jump(config::Config, polytopes::Vector{Polyhedron{T}}) where T
    
    polytopes_jump = Vector{Model()}()
    for polytope in polytopes
        polytope_jump = polyhedra_to_jump(config, polytope)
        push!(polytopes_jump, polytope_jump)
    end
    return polytopes_jump
end

"""
    check_intersection(intersecting_polytopes)

Check that the intersection of `intersecting_polytopes` is non-empty (prints the intersection size).

Note: this path has caused numerical issues with `Polyhedra.jl`/`CDDLib` in the past; use with care.
"""
function check_intersection(intersecting_polytopes::Vector{Polyhedron})
    
    # Get the number of points in the intersection of all polytopes
    intersection_size = npoints(intersect(intersecting_polytopes...))
    
    # Print intersection size for debugging
    println("Intersection of all polytopes is a polytope with $intersection_size vertices")
    
    # Assert that the intersection is not empty
    @assert intersection_size ≥ 1 "There must be at least one point in the intersection of all polytopes"
end

"""
    compute_distance(config, lmo_list) -> (last_iterate, last_lmo_vertex, primal, fw_gap)

Multiple dispatch: compute the distance objective between `k` polytopes by running FW with
`config.max_iterations_opt`.
"""
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
"""
    compute_distance(config, lmo_list) -> (last_iterate, last_lmo_vertex, primal, fw_gap)

Multiple dispatch for `Vector{FrankWolfe.ConvexHullLMO}`.
"""
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
"""
    compute_distance(config, lmo_list) -> (last_iterate, last_lmo_vertex, primal, fw_gap)

Multiple dispatch for `Vector{FrankWolfe.MathOptLMO}`.
"""
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

"""
    save_polytopes(filename, vertices, shifted_vertices, primal, fw_gap)

Save instance data to a `.jld2` file.
"""
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

"""
    save_polytopes(config, vertices, shifted_vertices, primal, fw_gap) -> filename

Multiple dispatch: save instance data to an auto-generated `.jld2` filename.
"""
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

"""
    load_polytopes(filename) -> (vertices, shifted_vertices, primal, fw_gap)

Load instance data from a `.jld2` file.
"""
function load_polytopes(filename::String)
    
    f = load(filename)
    vertices = f["vertices"]
    shifted_vertices = f["shifted_vertices"]
    primal = f["primal"]
    fw_gap = f["fw_gap"]
    
    return vertices, shifted_vertices, primal, fw_gap
end
