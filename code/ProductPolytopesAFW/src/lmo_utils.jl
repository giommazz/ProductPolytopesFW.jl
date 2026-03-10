# `lmos.jl`
"""
    MatrixConvexHullLMO(vertices; cache_cap=nothing)

Memory-efficient convex-hull LMO backed by a vertex matrix.

Vertices are stored row-wise in a matrix of size (`nverts x dim`), each row is one candidate atom.
Using a matrix representation is advantageous to use BLAS operations and handle Vector{Vector}.

Given a direction and a vertex matrix, LMO searches over rows/vertices and returns dense Vector
while reusing cache and buffers to reduce reallocations and time. 

The (capped) vertex cache stores dense copies of vertices recently selected by the LMO, so repeated queries 
can reuse them instead of copying the same row from the matrix again (reduce reallocations).
Clock (second-chance) policy → when cache is full, cycle over cache slots:
- if a slot was not used recently, it is reused immediately;
- if it was used recently, it gets one extra chance and is skipped once.

Note that, on convex hulls with dense vertex matrices and a lot of vertices,
MatrixConvexHullLMO with use_optimized_search=true should be systematically faster 
than ConvexHullLMO (which is generic by design)

"""
mutable struct MatrixConvexHullLMO{T,MT<:AbstractMatrix{T}} <: FrankWolfe.LinearMinimizationOracle
    vertices::MT                    # Row-wise vertex table (nverts x dim)
    cache_cap::Int                  # Max number of cached materialized vertices
    use_optimized_search::Bool      # true: efficient GEMV-based score computation, false: row-by-row scan
    idx_to_slot::Vector{Int}        # Vertex index -> cache slot (0 means "not cached")
    slot_to_idx::Vector{Int}        # Cache slot -> vertex index (0 means "empty")
    slot_refbit::BitVector          # true: recently referenced/used, false: first candidate to be replaced
    slot_store::Vector{Vector{T}}   # Cached dense vertices
    cache_len::Int                  # Number of currently filled slots
    clock_hand::Int                 # Current pointer in clock replacement scan
    score_work::Vector{T}           # Temporary buffer, scratch scores for one matrix-vector product (length nverts)
    direction_work::Vector{T}       # Temporary buffer, scratch contiguous direction for BLAS-friendly mul!
end

const MATRIX_LMO_CACHE_TARGET_BYTES = 64 * 1024 * 1024 # 64 MiB budget per LMO cache

"""
    CopyExtremePointLMO(lmo)

Wrapper LMO that "materializes" the returned extreme point as a dense `Vector`.

Why:
- `FrankWolfe.ConvexHullLMO` is often built from a list of vertex "views" (e.g. `eachrow(V)`),
  so `compute_extreme_point` returns a `SubArray`.
- If using `FrankWolfe.BlockVector{Float64, Vector{Float64}, ...}` iterates
  (see `find_starting_point`), the active-set expects atoms with dense `Vector` blocks.

This wrapper keeps the memory benefit of storing vertices as views, while ensuring
type-stable atoms compatible with the active-set machinery of ConvexHullLMO, by copying 
the selected vertex on output.
"""
struct CopyExtremePointLMO{L} <: FrankWolfe.LinearMinimizationOracle
    lmo::L
end

function FrankWolfe.compute_extreme_point(
    lmo::CopyExtremePointLMO,
    direction;
    v=nothing,
    kwargs...,
)
    s = FrankWolfe.compute_extreme_point(lmo.lmo, direction; v=v, kwargs...)
    return Vector(s)
end

"""
    _matrix_lmo_default_cache_cap(nverts, dim, T)

Heuristic default cap for cached materialized vertices, based on a fixed byte
budget per LMO. The cap is then clamped to `[32, nverts]`.
"""
function _matrix_lmo_default_cache_cap(nverts::Int, dim::Int, ::Type{T}) where {T}
    # Approximate bytes needed to store one dense vertex.
    bytes_per_vertex = max(dim * sizeof(T), 1)
    # Fit as many vertices as possible into the configured budget.
    # `÷` is integer division (quotient), so the cap is an integer count of vertices.
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
skip it once, pick the first slot with `refbit=false`.

`refbit` means "recently referenced": true if the slot was touched recently,
false if it has not been used since the previous clock pass.
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

Return cached dense vertex for `idx` if present, otherwise `nothing`.
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

Insert/update cached materialized `vertex` for row index `idx` into cache.
If the cache is full, reuse one older slot chosen by the "clock rule":
scan cache slots in a circular order and take the first one that was not
referenced recently.
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

Uses one matrix-vector product with preallocated buffers.
Return the row index minimizing `<V[i,:], direction>` for `V = lmo.vertices`.
"""
@inline function _best_row_index_optimized!(lmo::MatrixConvexHullLMO{T}, direction) where {T}
    # Copy `direction` into preallocated contiguous temp (reusable) buffer
    # → `mul!` can use fast path and avoid fresh allocations.
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

Row-by-row search: scan rows and evaluate one dot product per row.
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

    # Find the index of the best vertex in the atom matrix.
    best_idx = lmo.use_optimized_search ?
        # Score all rows efficiently at once via matrix-vector multiplication.
        _best_row_index_optimized!(lmo, direction) :
        # Inspect rows one by one with explicit dot products (less efficient)
        _best_row_index_rowscan(lmo, direction)

    # Try to reuse a previously materialized dense vertex.
    cached_vertex = _cache_lookup(lmo, best_idx)
    if v !== nothing && length(v) != dim
        throw(DimensionMismatch("Output buffer has length $(length(v)), expected $dim."))
    end

    # Return dense vectors so all atoms share the same block type in AFW active sets.
    if cached_vertex === nothing # vertex is not stored in cache yet
        # No output buffer (to write new extreme point) given: materialize once and cache it.
        if v === nothing
            new_vertex = Vector{T}(undef, dim)
            _copy_row!(new_vertex, V, best_idx)
            _cache_insert!(lmo, best_idx, new_vertex)
            return new_vertex
        end
        # Fill given buffer (`v != nothing`), and cache a dedicated copy.
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
    create_product_lmo(lmo_list::AbstractVector{L}) where {L<:FrankWolfe.LinearMinimizationOracle}

Build and return a `FrankWolfe.ProductLMO` from a list of LMOs.
"""
function create_product_lmo(lmo_list::AbstractVector{L}) where {L<:FrankWolfe.LinearMinimizationOracle}
    lmos_tuple = Tuple(lmo_list)
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
If `config.cvxhflag` is true, the convex-hull backend is selected by
`config.convex_hull_backend`:
- `"matrix"`: `MatrixConvexHullLMO`
- `"vector"`: `FrankWolfe.ConvexHullLMO` built from row views (best effort)
If `config.cvxhflag` is false, returns `FrankWolfe.MathOptLMO` objects
built from JuMP backends.
"""
# Initialize LMOs for a list of polytope vertex matrices.
function create_lmos(config::Config, vertices::Vector{Matrix{T}}) where T
    # Initialize data structures
    if config.cvxhflag lmo_list = Vector{FrankWolfe.LinearMinimizationOracle}() else lmo_list = Vector{FrankWolfe.MathOptLMO}() end

    # Create LMOs
    for V in vertices
        if config.cvxhflag
            if config.convex_hull_backend == "matrix"
                cache_cap = config.matrix_lmo_cache_cap == -1 ? nothing : config.matrix_lmo_cache_cap
                lmo = MatrixConvexHullLMO(
                    V;
                    cache_cap=cache_cap,
                    use_optimized_search=config.matrix_lmo_use_optimized_search,
                )
                push!(lmo_list, lmo)
            elseif config.convex_hull_backend == "vector"
                # Best effort for ConvexHullLMO:
                # - store vertex *views* to avoid duplicating `V`;
                # - but materialize each selected extreme point as a dense `Vector`,
                #   so AFW's active set sees consistent atom types.
                base_lmo = FrankWolfe.ConvexHullLMO(collect(eachrow(V)))
                lmo = CopyExtremePointLMO(base_lmo)
                push!(lmo_list, lmo)
            else
                error("Invalid convex_hull_backend='$(config.convex_hull_backend)'. Expected 'matrix' or 'vector'.")
            end
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
