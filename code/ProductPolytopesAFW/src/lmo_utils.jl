# `lmos.jl`
"""
    MatrixConvexHullLMO(vertices; cache_cap=nothing)

Memory-efficient convex-hull LMO backed by a vertex matrix.
It scans rows of `vertices` without materializing `Vector{Vector}`.

To reduce repeated allocations in AFW, it keeps a (capped) cache of materialized
vertices selected by the oracle. The cache uses a clock (second-chance) policy.
"""
mutable struct MatrixConvexHullLMO{T,MT<:AbstractMatrix{T}} <: FrankWolfe.LinearMinimizationOracle
    vertices::MT                    # Row-wise vertex table (nverts x dim)
    cache_cap::Int                  # Max number of cached materialized vertices
    use_optimized_search::Bool      # true: GEMV-based score computation, false: row-by-row scan
    idx_to_slot::Vector{Int}        # Vertex index -> cache slot (0 means "not cached")
    slot_to_idx::Vector{Int}        # Cache slot -> vertex index (0 means "empty")
    slot_refbit::BitVector          # true: recently used, false: first candidate to be replaced
    slot_store::Vector{Vector{T}}   # Cached dense vertices
    cache_len::Int                  # Number of currently filled slots
    clock_hand::Int                 # Current pointer in clock replacement scan
    score_work::Vector{T}           # Scratch scores for one matrix-vector product (length nverts)
    direction_work::Vector{T}       # Scratch contiguous direction for BLAS-friendly mul!
end

const MATRIX_LMO_CACHE_TARGET_BYTES = 64 * 1024 * 1024 # 64 MiB budget per LMO cache

"""
    _matrix_lmo_default_cache_cap(nverts, dim, T)

Heuristic default cap for cached materialized vertices, based on a fixed byte
budget per LMO. The cap is then clamped to `[32, nverts]`.
"""
function _matrix_lmo_default_cache_cap(nverts::Int, dim::Int, ::Type{T}) where {T}
    # Approximate bytes needed to store one dense vertex.
    bytes_per_vertex = max(dim * sizeof(T), 1)
    # Fit as many vertices as possible into the configured budget.
    raw_cap = MATRIX_LMO_CACHE_TARGET_BYTES ÷ bytes_per_vertex
    return clamp(raw_cap, 32, nverts)
end

"""
    MatrixConvexHullLMO(vertices; cache_cap=nothing, use_optimized_search=true)

Build a matrix-backed convex-hull LMO with a capped vertex cache.

- `cache_cap = nothing`: use an automatic cap from `_matrix_lmo_default_cache_cap`.
- `cache_cap = 0`: disable caching.
- `cache_cap > 0`: keep at most that many materialized vertices in cache.
- `use_optimized_search = true`: pick extreme points via GEMV (`mul!`) + argmin.
- `use_optimized_search = false`: use row-by-row dot scan (legacy path).
"""
function MatrixConvexHullLMO(
    vertices::MT;
    cache_cap::Union{Nothing,Int}=nothing,
    use_optimized_search::Bool=true,
) where {T,MT<:AbstractMatrix{T}}
    nverts, dim = size(vertices)
    if nverts == 0
        error("MatrixConvexHullLMO received an empty vertex matrix.")
    end
    if cache_cap !== nothing && cache_cap < 0
        error("`cache_cap` must be >= 0, got $cache_cap.")
    end
    # If unspecified, compute a cap from byte budget and vertex size.
    cap = cache_cap === nothing ? _matrix_lmo_default_cache_cap(nverts, dim, T) : cache_cap
    cap = clamp(cap, 0, nverts)
    return MatrixConvexHullLMO{T,MT}(
        vertices,
        cap,
        use_optimized_search,
        zeros(Int, nverts),
        zeros(Int, cap),
        falses(cap),
        Vector{Vector{T}}(undef, cap),
        0,
        1,
        zeros(T, nverts),
        zeros(T, dim),
    )
end

"""
    _clock_pick_slot_to_replace!(lmo)

Choose one cache slot to reuse when the cache is full.
Clock rule: if a slot was used recently (`refbit=true`), clear its flag and
skip it once; pick the first slot with `refbit=false`.
"""
@inline function _clock_pick_slot_to_replace!(lmo::MatrixConvexHullLMO)
    # Second-chance policy:
    # - refbit=false => choose this slot now
    # - refbit=true  => clear flag and give one extra turn
    while true
        slot = lmo.clock_hand
        if !lmo.slot_refbit[slot]
            lmo.clock_hand = slot == lmo.cache_cap ? 1 : slot + 1
            return slot
        end
        lmo.slot_refbit[slot] = false
        lmo.clock_hand = slot == lmo.cache_cap ? 1 : slot + 1
    end
end

"""
    _cache_lookup(lmo, idx)

Return cached dense vertex for `idx` if present; otherwise `nothing`.
If `idx` is already in cache, mark its slot as recently used.
"""
@inline function _cache_lookup(lmo::MatrixConvexHullLMO{T}, idx::Int) where {T}
    slot = lmo.idx_to_slot[idx]
    if slot == 0
        return nothing
    end
    # Mark as recently used so this slot is less likely to be replaced soon.
    lmo.slot_refbit[slot] = true
    return lmo.slot_store[slot]::Vector{T}
end

"""
    _cache_insert!(lmo, idx, vertex)

Insert/update cached materialized `vertex` for row index `idx`.
If the cache is full, reuse one older slot chosen by the clock rule.
"""
@inline function _cache_insert!(lmo::MatrixConvexHullLMO{T}, idx::Int, vertex::Vector{T}) where {T}
    # Fast path: cache disabled.
    if lmo.cache_cap == 0
        return
    end
    slot = if lmo.cache_len < lmo.cache_cap
        # Still room: append into first unused slot.
        lmo.cache_len += 1
        lmo.cache_len
    else
        # Full cache: reuse one old slot with clock policy.
        slot = _clock_pick_slot_to_replace!(lmo)
        old_idx = lmo.slot_to_idx[slot]
        if old_idx != 0
            # Clear reverse mapping for the replaced vertex.
            lmo.idx_to_slot[old_idx] = 0
        end
        slot
    end
    lmo.slot_store[slot] = vertex
    lmo.slot_to_idx[slot] = idx
    lmo.slot_refbit[slot] = true
    lmo.idx_to_slot[idx] = slot
end

"""
    _copy_row!(dest, V, row)

Copy row `row` of matrix `V` into dense vector `dest`.
"""
@inline function _copy_row!(dest::AbstractVector{T}, V::AbstractMatrix{T}, row::Int) where {T}
    @inbounds for j in eachindex(dest)
        dest[j] = V[row, j]
    end
    return dest
end

"""
    _best_row_index_optimized!(lmo, direction)

Return the row index minimizing `<V[i,:], direction>` for `V = lmo.vertices`.
Uses one matrix-vector product with preallocated scratch buffers.
"""
@inline function _best_row_index_optimized!(lmo::MatrixConvexHullLMO{T}, direction) where {T}
    # Copy into contiguous scratch so `mul!` can use a fast path.
    @inbounds for j in eachindex(lmo.direction_work)
        lmo.direction_work[j] = direction[j]
    end

    # Compute all row scores in one pass: score[i] = <V[i,:], direction>.
    LinearAlgebra.mul!(lmo.score_work, lmo.vertices, lmo.direction_work)

    # Argmin on the precomputed scores.
    best_idx = 1
    best_val = lmo.score_work[1]
    @inbounds for i in 2:length(lmo.score_work)
        val = lmo.score_work[i]
        if val < best_val
            best_val = val
            best_idx = i
        end
    end
    return best_idx
end

"""
    _best_row_index_rowscan(lmo, direction)

Legacy row-by-row search: scan rows and evaluate one dot product per row.
"""
@inline function _best_row_index_rowscan(lmo::MatrixConvexHullLMO, direction)
    V = lmo.vertices
    best_idx = 1
    best_val = dot(view(V, 1, :), direction)
    @inbounds for i in 2:size(V, 1)
        val = dot(view(V, i, :), direction)
        if val < best_val
            best_val = val
            best_idx = i
        end
    end
    return best_idx
end

"""
    FrankWolfe.compute_extreme_point(lmo::MatrixConvexHullLMO, direction; v=nothing, kwargs...)

Return the vertex `s` in `lmo.vertices` that minimizes `<s, direction>`.

- Use one matrix-vector product to score all rows and pick the best index.
- Reuse a cached dense vertex when available.
- If `v === nothing`, returns a dense vector (possibly the cached object).
- If `v` is provided, writes into `v` and returns `v`.

Cache keeps object identity stable for repeated vertex selections, helpful for AFW active-set management.
"""
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

    # Oracle step: choose optimized GEMV path or legacy rowscan path.
    best_idx = lmo.use_optimized_search ?
        _best_row_index_optimized!(lmo, direction) :
        _best_row_index_rowscan(lmo, direction)

    # Try to reuse a previously materialized dense vertex.
    cached_vertex = _cache_lookup(lmo, best_idx)
    if v !== nothing && length(v) != dim
        throw(DimensionMismatch("Output buffer has length $(length(v)), expected $dim."))
    end

    # Return dense vectors so all atoms share the same block type in AFW active sets.
    if cached_vertex === nothing
        if v === nothing
            # No output buffer: materialize once and cache it.
            new_vertex = Vector{T}(undef, dim)
            _copy_row!(new_vertex, V, best_idx)
            _cache_insert!(lmo, best_idx, new_vertex)
            return new_vertex
        end
        # Buffer path: fill caller buffer, and cache a dedicated copy.
        _copy_row!(v, V, best_idx)
        _cache_insert!(lmo, best_idx, copy(v))
        return v
    end

    if v === nothing
        # Cache hit: return the stable cached object (important for AFW atom identity).
        return cached_vertex
    end
    # Cache hit with output buffer requested.
    copyto!(v, cached_vertex)
    return v
end

"""
    create_product_lmo(lmo_list::Vector{FrankWolfe.LinearMinimizationOracle})

Build and return a `FrankWolfe.ProductLMO` from a list of generic LMOs.
"""
function create_product_lmo(lmo_list::Vector{FrankWolfe.LinearMinimizationOracle})
    
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmo_list)
    
    # Create and return a ProductLMO object
    return FrankWolfe.ProductLMO(lmos_tuple)
end

"""
    create_product_lmo(lmo_list::Vector{FrankWolfe.ConvexHullLMO})

Build and return a `FrankWolfe.ProductLMO` from a list of `ConvexHullLMO`s.
"""
# (Multiple dispatch)
function create_product_lmo(lmo_list::Vector{FrankWolfe.ConvexHullLMO})
    
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmo_list)
    
    # Create and return a ProductLMO object
    return FrankWolfe.ProductLMO(lmos_tuple)
end

"""
    create_product_lmo(lmo_list::Vector{FrankWolfe.MathOptLMO})

Build and return a `FrankWolfe.ProductLMO` from a list of `MathOptLMO`s.
"""
# (Multiple dispatch)
function create_product_lmo(lmo_list::Vector{FrankWolfe.MathOptLMO})
    
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmo_list)
    
    # Create and return a ProductLMO object
    return FrankWolfe.ProductLMO(lmos_tuple)
end

"""
    find_starting_point(config, prod_lmo)

Create an initial `FrankWolfe.BlockVector` by querying one extreme point per
block LMO and materializing each block as a dense `Vector{Float64}`.
"""
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

"""
    create_lmos(config, vertices::Vector{Matrix{T}})

Create one LMO per polytope matrix in `vertices`.
If `config.cvxhflag` is true, returns matrix-backed convex-hull LMOs;
otherwise returns `FrankWolfe.MathOptLMO` objects built from JuMP backends.
"""
# Initialize LMOs for given sets of `vertices` k=1 polytope
# Depending on `cvxhflag`, create either `FrankWolfe.ConvexHullLMO` (true) or `FrankWolfe.MathOptLMO` (false) objects.
function create_lmos(config::Config, vertices::Vector{Matrix{T}}) where T
    # Initialize data structures
    if config.cvxhflag lmo_list = Vector{FrankWolfe.LinearMinimizationOracle}() else lmo_list = Vector{FrankWolfe.MathOptLMO}() end

    # Create LMOs
    for V in vertices
        if config.cvxhflag # Matrix-backed convex hull LMOs (memory-safe and AFW-compatible)
            cache_cap = config.matrix_lmo_cache_cap == -1 ? nothing : config.matrix_lmo_cache_cap
            lmo = MatrixConvexHullLMO(
                V;
                cache_cap=cache_cap,
                use_optimized_search=config.matrix_lmo_use_optimized_search,
            )
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

"""
    create_lmos(config, vertices::Vector{Vector{Matrix{T}}})

Create LMOs for multiple polytope collections, returning one LMO list per
entry of `vertices`.
"""
# (Multiple dispatch) k polytopes
function create_lmos(config::Config, vertices::Vector{Vector{Matrix{T}}}) where T
    lmo_list = []
    for vs in vertices
        lmos = create_lmos(config, vs)
        push!(lmo_list, lmos)
    end
    return lmo_list
end
