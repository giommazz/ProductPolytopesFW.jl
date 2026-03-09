# `polytope_generation.jl`

"""
    compute_anchor(config, vertices)

Compute the anchor point in `R^n` according to `config.intersection_anchor`.
"""
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

"""
    generate_nonintersecting_bounds(config; margin=10.0, stepsize=100.0, start_point=-100.0)

Generate `k` axis-aligned boxes in `R^n` with guaranteed separation (non-overlap).
"""
function generate_nonintersecting_bounds(config::Config; margin::Float64=10.0, stepsize::Float64=100.0, start_point::Float64=-100.0)
    
    bounds_list = Vector{Vector{Tuple{Float64, Float64}}}(undef, config.k)

    # Generate bounds for P₁ randomly as pairs (LB, UB), within given range
    lower_bounds = [start_point + stepsize * rand(Float64) for _ in 1:config.n]
    upper_bounds = [lower_bounds[d] + (stepsize - lower_bounds[d]) * rand(Float64) for d in 1:config.n]
    bounds_list[1] = [(lower_bounds[d], upper_bounds[d]) for d in 1:config.n]

    # Generate bounds for P₂, ..., Pₖ while ensuring no overlap
    for i in 2:config.k
        # Current polytope
        # Note that,at start of iteration `i`, `upper_bounds` contains UBs of previous polytope Pᵢ₋₁
        lower_bounds = [upper_bounds[d] + margin + stepsize * rand(Float64) for d in 1:config.n]
        upper_bounds = [lower_bounds[d] + stepsize * rand(Float64) for d in 1:config.n]
        bounds_list[i] = [(lower_bounds[d], upper_bounds[d]) for d in 1:config.n]
    end
    
    return bounds_list
end

"""
    generate_polytope(config, idx, bounds) -> (vertices, p1_vertex_mean_or_nothing)

Generate one polytope by sampling `config.n_points[idx]` points within the coordinate-wise `bounds`,
using `config.vertex_sampling`.
"""
function generate_polytope(config::Config, idx::Int, bounds::Vector{Tuple{T, T}}) where T
    
    # Check that the number of columns specified in `config.n` matches the length of `bounds`
    if config.n != length(bounds)
        error("Mismatch: The number of columns in `config` (config.n = $(config.n)) must match the number of elements in `bounds` (length = $(length(bounds))).")
    end

    # Initialize matrix of zeros with given number of points (rows) and dimensions (cols)
    vertices = zeros(Float64, config.n_points[idx], config.n)

    if config.vertex_sampling == "box_uniform"
        # "box_uniform": sample points uniformly inside the axis-aligned box.
        # -> can be anisotropic across coordinates because each dimension has its own (LB, UB) range.
        for i in 1:config.n_points[idx]
            vertices[i, :] = [bounds[d][1] + (bounds[d][2] - bounds[d][1]) * rand(Float64) for d in 1:config.n]
        end
    elseif config.vertex_sampling == "sphere"
        # "sphere": sample points on a sphere fully and strictly contained in the box.
        # -> produces more isotropic point clouds while preserving non-intersection:
        # -> convex hull stays fully and strictly inside its box, boxes are disjoint by construction.
        center = [(bounds[d][1] + bounds[d][2]) / 2 for d in 1:config.n]
        box_half_lengths = [(bounds[d][2] - bounds[d][1]) / 2 for d in 1:config.n]
        radius = config.sphere_radius_factor * minimum(box_half_lengths)
        for i in 1:config.n_points[idx]
            g = randn(Float64, config.n) # random direction (Gaussian N(0,1)) for current point
            ng = norm(g)
            while ng == 0.0 # generate nonzero points
                g = randn(Float64, config.n)
                ng = norm(g)
            end
            # each point is guaranteed to be *strictly* inside the sphere: define the following quantities
            # u = g/ng of unit norm;        R_i = ρ * min_d hᵢ[d];          v = cᵢ + Rᵢ * u,
            # then
            # |v[d] - cᵢ[d]| = Rᵢ * |u[d]| ≤ Rᵢ ≤ ρ * hᵢ[d] < hᵢ[d]
            vertices[i, :] = center .+ (radius / ng) .* g # point on sphere of radius `radius`
        end
    elseif config.vertex_sampling == "ellipsoid"
        # "ellipsoid": sample points on an axis-aligned ellipsoid fully and strictly contained in the box.
        # -> avoids the "min side-length" collapse of the sphere construction while staying inside the box:
        # -> v = c + ρ * (h .* u), where h are the box half-lengths and u is a random unit direction.
        center = [(bounds[d][1] + bounds[d][2]) / 2 for d in 1:config.n]
        box_half_lengths = [(bounds[d][2] - bounds[d][1]) / 2 for d in 1:config.n]
        axes = config.sphere_radius_factor .* box_half_lengths
        for i in 1:config.n_points[idx]
            g = randn(Float64, config.n) # random direction (Gaussian N(0,1)) for current point
            ng = norm(g)
            while ng == 0.0 # generate nonzero points
                g = randn(Float64, config.n)
                ng = norm(g)
            end
            u = g ./ ng # random unit direction
            vertices[i, :] = center .+ axes .* u
        end
    else
        error("Unknown vertex sampling strategy: $(config.vertex_sampling)")
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

"""
    generate_nonintersecting_polytopes(config, bounds) -> (vertices, p1_vertex_mean, primal, fw_gap)

Generate `k` non-intersecting polytopes inside the provided `bounds`, then certify non-intersection by
computing the distance objective via FW.
"""
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

"""
    intersect_polytopes(config, vertices) -> shifted_vertices

Shift the input polytopes so that their intersection contains the anchor point specified by `config`.
"""
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

"""
    generate_polytopes(config) -> (vertices, shifted_vertices, primal, fw_gap)

Generate `k` non-intersecting polytopes, then shift them according to the anchor/reference-point rules in
`config` so that the shifted polytopes intersect.
"""
function generate_polytopes(config::Config)

    # `n_points` contains n. of vertices used to generate each polytope
    if config.verbose
        println("Generating $(config.k) intersecting polytopes of dimension $(config.n) (n. vertices for each: $(config.n_points))")
    end
    
    # Generate random, non-intersecting boxes, then a polytope in each of the boxes
    bounds_list = generate_nonintersecting_bounds(config)

    if config.verbose
        fmt(x::Float64) = isfinite(x) ? round(x; digits=6) : x

        if config.vertex_sampling == "sphere"
            # Diagnostics: quantify how small the inscribed sphere radius is relative to the box size.
            # For each polytope box i, let hᵢ[d] = (UBᵢ[d]-LBᵢ[d])/2 be the half-lengths, rᵢ = ρ*min_d hᵢ[d],
            # and define shrinkᵢ = rᵢ / ‖hᵢ‖ (fraction of the box half-diagonal covered by the sphere radius).
            rho = config.sphere_radius_factor
            println("\nSphere sampling box stats (ρ=$(fmt(rho))) (rounded to 6 decimals):")
            println("  Legend / definitions:")
            println("    h[d]     = (UB[d]-LB[d])/2                         (box half-lengths)")
            println("    h_min    = min_d h[d], h_med = median_d h[d], h_mean = mean_d h[d]")
            println("    min/med  = h_min / h_med, max/min = max_d h[d] / h_min")
            println("    r        = ρ * h_min                               (sphere radius)")
            println("    |h|      = ‖h‖₂ = sqrt(∑_d h[d]^2)                 (box half-diagonal)")
            println("    shrink   = r / |h|                                 (radius relative to box half-diagonal)")
            println("  Values:")
            shrinks = Float64[]
            for (i, b) in enumerate(bounds_list)
                h = [(ub - lb) / 2 for (lb, ub) in b]
                hmin = minimum(h)
                hmed = median(h)
                hmean = mean(h)
                hmax = maximum(h)
                min_over_med = hmed == 0.0 ? NaN : hmin / hmed
                max_over_min = hmin == 0.0 ? Inf : hmax / hmin
                r = rho * hmin
                hd = norm(h)
                shrink = hd == 0.0 ? NaN : r / hd
                push!(shrinks, shrink)
                println("  P$i: h_min=$(fmt(hmin)) h_med=$(fmt(hmed)) h_mean=$(fmt(hmean)) min/med=$(fmt(min_over_med)) max/min=$(fmt(max_over_min)) r=$(fmt(r)) shrink=$(fmt(shrink))")
            end
            println("  shrink across polytopes: min=$(fmt(minimum(shrinks))) mean=$(fmt(mean(shrinks))) max=$(fmt(maximum(shrinks)))\n")
        elseif config.vertex_sampling == "ellipsoid"
            # Diagnostics: quantify how the inscribed ellipsoid compares to the box size.
            # For each polytope box i, let hᵢ[d] = (UBᵢ[d]-LBᵢ[d])/2 be the half-lengths, and aᵢ[d] = ρ*hᵢ[d] the ellipsoid semi-axes.
            # We report a_min = min_d aᵢ[d] and shrink_min = a_min / ‖hᵢ‖, plus shrink_diag = ‖aᵢ‖ / ‖hᵢ‖ (= ρ).
            rho = config.sphere_radius_factor
            println("\nEllipsoid sampling box stats (ρ=$(fmt(rho))) (rounded to 6 decimals):")
            println("  Legend / definitions:")
            println("    h[d]        = (UB[d]-LB[d])/2                       (box half-lengths)")
            println("    a[d]        = ρ * h[d]                              (ellipsoid semi-axis lengths)")
            println("    h_min/h_med/h_mean and min/med/max/min as for sphere")
            println("    a_min       = min_d a[d] = ρ * h_min                (smallest semi-axis)")
            println("    |h|         = ‖h‖₂ = sqrt(∑_d h[d]^2)               (box half-diagonal)")
            println("    shrink_min  = a_min / |h|                           (smallest semi-axis relative to box half-diagonal)")
            println("    shrink_diag = ‖a‖₂ / ‖h‖₂ = ρ                       (ellipsoid half-diagonal relative to box half-diagonal)")
            println("  Values:")
            shrinks_min = Float64[]
            for (i, b) in enumerate(bounds_list)
                h = [(ub - lb) / 2 for (lb, ub) in b]
                hmin = minimum(h)
                hmed = median(h)
                hmean = mean(h)
                hmax = maximum(h)
                min_over_med = hmed == 0.0 ? NaN : hmin / hmed
                max_over_min = hmin == 0.0 ? Inf : hmax / hmin

                a_min = rho * hmin
                hd = norm(h)
                shrink_min = hd == 0.0 ? NaN : a_min / hd
                shrink_diag = rho
                push!(shrinks_min, shrink_min)

                println("  P$i: h_min=$(fmt(hmin)) h_med=$(fmt(hmed)) h_mean=$(fmt(hmean)) min/med=$(fmt(min_over_med)) max/min=$(fmt(max_over_min)) a_min=$(fmt(a_min)) shrink_min=$(fmt(shrink_min)) shrink_diag=$(fmt(shrink_diag))")
            end
            println("  shrink_min across polytopes: min=$(fmt(minimum(shrinks_min))) mean=$(fmt(mean(shrinks_min))) max=$(fmt(maximum(shrinks_min)))\n")
        end
    end
    vertices, p1_vertex_mean, primal, fw_gap = generate_nonintersecting_polytopes(config, bounds_list)

    # Intersect them according to `intersection.anchor` and `intersection.reference_point` in `config`
    shifted_vertices = intersect_polytopes(config, vertices)

    return vertices, shifted_vertices, primal, fw_gap
end
