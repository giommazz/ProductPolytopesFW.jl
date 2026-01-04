# `polytope_generation.jl`

# Compute the anchor point according to `config.intersection_anchor`
function compute_anchor(config::Config, vertices::Vector{Matrix{T}}) where T
    if isempty(vertices)
        error("Cannot compute anchor: no polytopes provided.")
    end

    anchor_type = config.intersection_anchor
    rng = MersenneTwister(config.seed) # `MersenneTwister`: pseudorand num. generator, seeded for reproducible randomness with `config.seed`
    Tvert = eltype(vertices[1]) # element type of vertex coordinates

    if anchor_type == "p1_center" # Vertex mean of the first polytope P₁
        return vertex_mean(vertices[1])
    elseif anchor_type == "p1_vertex" # Random vertex of P₁
        return random_vertex(vertices[1], rng)
    elseif anchor_type == "p1_random" # Random convex combination of vertices of P₁
        return random_convex_combination(vertices[1], rng)
    elseif anchor_type == "origin" # Zero vector in Rⁿ
        return zeros(Tvert, config.n)
    elseif anchor_type == "global_mean" # Mean of per-polytope vertex means across all polytopes
        return global_mean_of_centers(config, vertices)
    elseif anchor_type == "random" # Random point in small interval [-1, 1]ⁿ (keep translations bounded + avoid huge coordinates)
        return 2 .* rand(rng, Tvert, config.n) .- one(Tvert)
    else
        error("Unknown intersection anchor: $(anchor_type)")
    end
end

# Function to generate non-overlapping bounds for multiple dimensions within [1e-03, 1e+03]
# `margin`: ensures min distance between current polytope's UB and next polytope's LB, to create non-overlapping polytopes
function generate_nonintersecting_bounds(config::Config; margin::Float64=10.0, stepsize::Float64=100.0, start_point::Float64=-100.0)
    
    bounds_list = Vector{Vector{Tuple{Float64, Float64}}}(undef, config.k)

    # Generate the first bounds randomly as pairs (LB, UB), within given range
    lower_bounds = [start_point + stepsize * rand(Float64) for _ in 1:config.n]
    upper_bounds = [lower_bounds[d] + (stepsize - lower_bounds[d]) * rand(Float64) for d in 1:config.n]
    bounds_list[1] = [(lower_bounds[d], upper_bounds[d]) for d in 1:config.n]

    # Generate subsequent bounds while ensuring no overlap
    for i in 2:config.k
        # Current polytope
        lower_bounds = [upper_bounds[d] + margin + stepsize * rand(Float64) for d in 1:config.n]
        upper_bounds = [lower_bounds[d] + stepsize * rand(Float64) for d in 1:config.n]
        bounds_list[i] = [(lower_bounds[d], upper_bounds[d]) for d in 1:config.n]
    end
    
    return bounds_list
end

# Generate one polytope: generate one random set of vertices, within given bounds for each dimension
function generate_polytope(config::Config, idx::Int, bounds::Vector{Tuple{T, T}}) where T
    
    # Check that the number of columns specified in `config.n` matches the length of `bounds`
    if config.n != length(bounds)
        error("Mismatch: The number of columns in `config` (config.n = $(config.n)) must match the number of elements in `bounds` (length = $(length(bounds))).")
    end

    # Initialize a matrix of zeros with the given number of points (rows) and columns
    vertices = zeros(Float64, config.n_points[idx], config.n)

    # Generate `config.n_points` vertices of size `config.n` within the bounds for each dimension
    for i in 1:config.n_points[idx]
        # Generate random Float64 within the given bounds for each coordinate
        vertices[i,:] = [bounds[d][1] + (bounds[d][2] - bounds[d][1]) * rand(Float64) for d in 1:config.n]
    end

    # Only compute vertex mean of P₁
    if idx == 1
        anc = vertex_mean(vertices)
    else
        anc = Nothing
    end

    # TODO: THIS CAUSES NUMERICAL ISSUES WITH THE Polyhedra.jl LIBRARY, DON'T USE FOR NOW (relevant functions commented)
    # Generate polytope as convex hull of given points, and remove redundant (i.e. non-vertex) points
    #vertices = nonredundant_polytope(vertices)
    
    return vertices, anc
end

# Generate k nonintersecting polytopes: generate k random sets of vertices, within given bounds for each dimension
function generate_nonintersecting_polytopes(config::Config, bounds::Vector{Vector{Tuple{T, T}}}) where T
    
    vertices = Vector{Matrix{T}}()
    p1_vertex_mean = Vector{T}()

    # Generate k polytopes within the specified bounds
    for i in 1:config.k
        if i == 1
            # Generate vertices and corresponding polytope
            verts, p1_vertex_mean = generate_polytope(config, i, bounds[i])
        else
            verts, _ = generate_polytope(config, i, bounds[i])
        end
        # Append to `vertices` and `vertex_means`
        push!(vertices, verts)
    end
    
    nonintersecting_polytopes_lmos = create_lmos(config, vertices)
    config_opt = modify_config(config, target_tolerance=config.target_tolerance_opt)    
    _, _, primal, fw_gap = compute_distance(config_opt, nonintersecting_polytopes_lmos)

    # Check that polytope intersection is empty
    if approxequal(primal, 0.0)
        error("Invalid polytopes: they should not intersect, but the total distance among them is $primal.")
    end

    return vertices, p1_vertex_mean, primal, fw_gap
end

# Function to intersect k polytopes so that their intersection contains the chosen anchor point
function intersect_polytopes(
    config::Config,
    vertices::Vector{Matrix{T}}
    ) where T
    
    if isempty(vertices)
        error("Cannot intersect polytopes: no polytopes provided.")
    end

    shifted_vertices = Vector{Matrix{T}}(undef, length(vertices))

    if config.intersection_anchor == "p1_p2_segment"
        # Legacy scheme: use closest pair (p₁, p₂) between P₁ and P₂, and vertex mean m₁ of P₁
        V1 = vertices[1]
        V2 = vertices[2]
        p1, p2 = closest_pair(config, V1, V2)
        m1 = vertex_mean(V1)
        t = config.intersection_anchor_t

        # Anchor lies on the segment p₁ → m₁ (or beyond if t > 1)
        anchor = p1 .+ t .* (m1 .- p1)

        # P₁: keep fixed
        shifted_vertices[1] = V1

        # P₂: reference point is p₂ from the closest pair
        shifted_vertices[2] = shift_polytope(V2, p2, anchor; stepsize=one(T))

        # P₃,…,P_k: reference points are vertices closest to the anchor
        for i in 3:length(vertices)
            ref_i = closest_pair(config, anchor, vertices[i])
            shifted_vertices[i] = shift_polytope(vertices[i], ref_i, anchor; stepsize=one(T))
        end
    else
        # Generic scheme: choose anchor, then compute reference points according to config.intersection_reference_point
        anchor = compute_anchor(config, vertices)
        refs = compute_reference_points(config, vertices, anchor)
        anchor_from_p1 = config.intersection_anchor in ["p1_center", "p1_vertex", "p1_random"]

        for i in 1:length(vertices)
            if i == 1 && anchor_from_p1
                # Anchor already lies in P₁; keep P₁ fixed
                shifted_vertices[i] = vertices[i]
            else
                shifted_vertices[i] = shift_polytope(vertices[i], refs[i], anchor; stepsize=one(T))
            end
        end
    end

    return shifted_vertices
end

# Generate nonintersecting polytopes, then intersect them according to the anchor/reference-point rules
function generate_polytopes(config::Config)

    # `n_points` contains n. of vertices used to generate each polytope
    println("Generating $(config.k) intersecting polytopes of dimension $(config.n) (n. vertices for each: $(config.n_points))")
    
    # Generate random, non-intersecting boxes, then a polytope in each of the boxes
    bounds_list = generate_nonintersecting_bounds(config)
    vertices, p1_vertex_mean, primal, fw_gap = generate_nonintersecting_polytopes(config, bounds_list)

    # Intersect them according to `intersection.anchor` and `intersection.reference_point` in `config`
    shifted_vertices = intersect_polytopes(config, vertices)

    return vertices, shifted_vertices, primal, fw_gap
end
