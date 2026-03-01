# `lmos.jl`
"""
    MatrixConvexHullLMO(vertices)

Memory-efficient convex-hull LMO backed by a vertex matrix.
It scans rows of `vertices` without materializing `Vector{Vector}`,
and returns a dense `Vector` extreme point (stable type for AFW active sets).
"""
struct MatrixConvexHullLMO{T,MT<:AbstractMatrix{T}} <: FrankWolfe.LinearMinimizationOracle
    vertices::MT
end

function FrankWolfe.compute_extreme_point(
    lmo::MatrixConvexHullLMO{T},
    direction;
    v=nothing,
    kwargs...,
) where T
    V = lmo.vertices
    nverts, dim = size(V)
    if nverts == 0
        error("MatrixConvexHullLMO received an empty vertex matrix.")
    end
    if length(direction) != dim
        throw(DimensionMismatch("Direction has length $(length(direction)), but vertices have dimension $dim."))
    end

    best_idx = 1
    best_val = dot(view(V, 1, :), direction)
    @inbounds for i in 2:nverts
        val = dot(view(V, i, :), direction)
        if val < best_val
            best_val = val
            best_idx = i
        end
    end

    # Return dense vectors so all atoms share the same block type in AFW active sets.
    if v === nothing
        return Vector{T}(view(V, best_idx, :))
    end
    if length(v) != dim
        throw(DimensionMismatch("Output buffer has length $(length(v)), expected $dim."))
    end
    @inbounds for j in 1:dim
        v[j] = V[best_idx, j]
    end
    return v
end

function create_product_lmo(lmo_list::Vector{FrankWolfe.LinearMinimizationOracle})
    
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmo_list)
    
    # Create and return a ProductLMO object
    return FrankWolfe.ProductLMO(lmos_tuple)
end
# (Multiple dispatch)
function create_product_lmo(lmo_list::Vector{FrankWolfe.ConvexHullLMO})
    
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmo_list)
    
    # Create and return a ProductLMO object
    return FrankWolfe.ProductLMO(lmos_tuple)
end
# (Multiple dispatch)
function create_product_lmo(lmo_list::Vector{FrankWolfe.MathOptLMO})
    
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmo_list)
    
    # Create and return a ProductLMO object
    return FrankWolfe.ProductLMO(lmos_tuple)
end

# Find starting point `x0` over the product of different LMOs
function find_starting_point(config::Config, prod_lmo::FrankWolfe.ProductLMO)
    # Prepare datafor `x0`, which is of type `FrankWolfe.BlockVector{Float64, Vector{Float64}, Tuple{Int64}}`
    # 1) Compute extreme points for each LMO in the product.
    #    Materialize as `Vector{Float64}` to avoid aliasing (e.g. SubArray views)
    #    since FrankWolfe algorithms may mutate `x0` in-place in `InplaceEmphasis()` mode.
    extreme_points = Vector{Vector{Float64}}(undef, length(prod_lmo.lmos))
    zero_dir = zeros(Float64, config.n)
    for (i, lmo) in enumerate(prod_lmo.lmos)
        s = FrankWolfe.compute_extreme_point(lmo, zero_dir)
        extreme_points[i] = Vector{Float64}(s) # copy/convert to dense Float64
    end
    # 2) Generate a vector of length `k`, where each entry is a tuple `(n,)`
    block_sizes = fill((config.n,), length(prod_lmo.lmos))
    # 3) Compute total size of `x0`
    total_size = sum([size[1] for size in block_sizes])

    # Instantiate `x0`
    return FrankWolfe.BlockVector(extreme_points, block_sizes, total_size)
end

# Initialize LMOs for given sets of `vertices` k=1 polytope
# Depending on `cvxhflag`, create either `FrankWolfe.ConvexHullLMO` (true) or `FrankWolfe.MathOptLMO` (false) objects.
function create_lmos(config::Config, vertices::Vector{Matrix{T}}) where T
    # Initialize data structures
    if config.cvxhflag lmo_list = Vector{FrankWolfe.LinearMinimizationOracle}() else lmo_list = Vector{FrankWolfe.MathOptLMO}() end

    # Create LMOs
    for V in vertices
        if config.cvxhflag # Matrix-backed convex hull LMOs (memory-safe and AFW-compatible)
            lmo = MatrixConvexHullLMO(V)
            # Update data structures
            push!(lmo_list, lmo)
        else    # FrankWolfe.MathOptLMO objects 
            # Instantiate Polyedra.Polyhedron objects
            poly = polytope(V)
            # Create JuMP.Model objects from Polyedra.Polyhedron objects
            jump_poly = polyhedra_to_jump(config, poly)
            # Ensure models are optimized
            optimize!(jump_poly)
            # Create FrankWolfe.MathOptLMO objects from JuMP.Model objects
            lmo = FrankWolfe.MathOptLMO(jump_poly.moi_backend)
            # Update data structures
            push!(lmo_list, lmo)
        end
    end
    return lmo_list
end
# (Multiple dispatch) k polytopes
function create_lmos(config::Config, vertices::Vector{Vector{Matrix{T}}}) where T
    lmo_list = []
    for vs in vertices
        lmos = create_lmos(config, vs)
        push!(lmo_list, lmos)
    end
    return lmo_list
end
