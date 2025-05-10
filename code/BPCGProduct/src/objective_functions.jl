# `objective_functions.jl`


# ------------------------------------------------------
# OBJECTIVE FUNCTION IMPLEMENTATIONS
# ------------------------------------------------------
"""
Function to compute form 1 of the objective function
f(x)  =     1/(2k) ∑_{1 ≤ i < j ≤ k} ‖xⁱ - xʲ‖           [form 1]
      =     1/(2k) [ (k-1) ∑ᵢ₌₁ᵏ‖xⁱ‖²  -  2 ∑_{1 ≤ i < j ≤ k} ⟨xⁱ, xʲ⟩ ]       [form 2.a]
      =     1/(2k) k∑ᵢ₌₁ᵏ‖xⁱ‖² - ‖∑_{1 ≤ i < j ≤ k} xⁱ‖²                    [form 2.b]
      =     1/(2k)<x, Mₖ x>, where Mₖ = (kI - 𝟏𝟏ᵀ), i.e., a (kn X kn) matrix with (k-1) on the diagonal and -1 elsewhere           [form 3]


ASYMPTOTIC COMPLEXITY OF THE FORMS
- forms 1 and 2.a have complexity O(k^2n) because they must compute Binomial(k, 2) = (k-1)k/2 terms, (one for each pair of different blocks), either
        squared Euclidean norms in ℝⁿ, or
        dot products in ℝⁿ
- form 2.b has complexity O(kn) because it only handles k blocks, each in ℝⁿ
- form 3 also has complexity O(kn) when we don't explicitly materialize the whole matrix. This can be achieved considering that
        Mₖx = kx - (𝟏ᵀx)𝟏 so we only need a linear number of calls

MEMORY OF THE FORMS (EXTRA RAM NEEDED BEYOND THE k BLOCKS)
- form 1: none
- form 2.b: one extra vector in ℝⁿ for the running sum
- form 3: same as 2.b if you apply Mₖ on the fly; huge (kn)² array if you store it
"""
# asymptotic complexity O(k²n), no extra memory complexity 
function convex_feasibility_objective_v1(x::FrankWolfe.BlockVector)
    
    sum_dist = 0.0
    # number of blocks
    k = length(x.block_sizes)
    for i in 1:k  
        for j in i+1:k
            xi = x.blocks[i]
            xj = x.blocks[j]
            curr = sum((xi - xj).^2)
            sum_dist += curr
        end
    end
    return 1.0 / (2k) * sum_dist
end

"""
Compute form 2
      f(x) = (1/(2k)) * [ k*Σ‖xᵢ‖²  -  ‖Σ xᵢ‖² ]
with a single pass through the `k` blocks.
"""
# asymptotic complexity O(kn), extra memory complexity O(n)
function convex_feasibility_objective_v2b(x::FrankWolfe.BlockVector)
    
    # number of blocks
    k = length(x.block_sizes)
    # dimension of each block
    n = length(x.blocks[1])
    # will contain Σ‖xᵢ‖²
    s_norm = 0.0
    # Allocate uninitialised vector in ℝⁿ to hold ∑xᵢ: `undef` is faster than `zeros` as it skips the default fill
    # basically: asks the Garbage Collector for n slots of *raw* memory and obtains them w/o touching the bytes
    # ⟹ no time spent overwriting whatever random bits were already in that RAM. So: allocation is light-fast
    s_vec = Vector{eltype(x)}(undef, n)
    # zero-initialize
    fill!(s_vec, 0)

    # fast, bounds‑check‑free loop
    # `@inbounds` (bounds check): every time Julia accesses v[i], it checks 1 ≤ i ≤ length(v) for segfaults
    #   `@inbounds` tells the compiler "I swear indices are safe, omit those checks"
    #   ⟹ machine code is shorter (no extra check) and may vectorize better
    # `@simd`: tells LLVM that no cross-iteration dependency exists ⟺ every iter can run independently
    @inbounds @simd for blk in x.blocks
        # add ‖⋅‖² of the current block to `s_norm`. 
        # `dot` is a BLAS Level‑1 call ⟺ highly optimised and parallel code on vector-vector operations
        s_norm += dot(blk, blk)
        # Performs the AXPY operation s_vec += 1.0⋅blk.
        # `axpy!` is another BLAS Level‑1 call ⟺ avoids temp allocations (≡ no GC overhead) and updates `s_vec` in-place
        LinearAlgebra.BLAS.axpy!(1.0, blk, s_vec)
    end
    # computes k⋅Σ‖xᵢ‖² - ‖∑xᵢ‖²
    return 0.5/k * (k*s_norm - dot(s_vec, s_vec))
end


# ------------------------------------------------------
# GRADIENT FUNCTION IMPLEMENTATIONS
# ------------------------------------------------------


# Compute the gradient (v1)
# ∇ⁱf(x) = 1/k [ (k-1)xⁱ - ∑_{j ≠ i}xʲ ]
# `storage` must be a FrankWolfe.BlockVector of the same shape as `x`
# asymptotic complexity O(k²n), extra memory complexity O(n)
function convex_feasibility_gradient_v1!(storage::FrankWolfe.BlockVector, x::FrankWolfe.BlockVector)
    
    k = length(x.block_sizes)
    n = Int(x.tot_size / k)
    for i in 1:k
        sum_terms = zeros(n)
        for j in 1:k
            if j != i
                sum_terms .+= x.blocks[j]
            end
        end
        storage.blocks[i] .= 1.0/k * ((k-1) * x.blocks[i] - sum_terms)
    end
end

# Compute the gradient (v2). Idea: compute s = ∑ᵏⱼ₌₁xʲ only once, because
# ∇ⁱf(x)    =   1/k [ (k-1)xⁱ - ∑_{j ≠ i}xʲ ]
#           =   1/k [ kxⁱ - xⁱ - ∑_{j ≠ i}xʲ ]
#           =   1/k [ kxⁱ - ∑ⱼ₌₁ᵏxʲ ]
#           =   1/k [ kxⁱ - s ]
# asymptotic complexity O(kn), extra memory complexity O(n)
function convex_feasibility_gradient_v2!(storage::FrankWolfe.BlockVector, x::FrankWolfe.BlockVector)
    k = length(x.block_sizes)
    n = length(x.blocks[1])
    # will contain ∑ⱼ₌₁ᵏ: one extra vector in ℝⁿ
    sum_vec = zeros(eltype(x), n)

    # build `sum_vec`
    @inbounds for blk in x.blocks
        # efficient implementation of sum_vec += blk
        BLAS.axpy!(one(eltype(x)), blk, sum_vec)
    end

    # reciprocal of k, i.e., 1/k
    scale = inv(k)
    # over blocks
    @inbounds @simd for i in eachindex(x.blocks)
        # over coordinates
        for j in 1:n
            storage.blocks[i][j] = scale * (k * x.blocks[i][j] - sum_vec[j])
        end
    end
end

# function used to generate a random FrankWolfe.BlockVector, to test convex_feasibility_objective variants
function random_blockvector(k::Integer, n::Integer; rng = Random.default_rng())
    blocks = [randn(rng, n) for _ in 1:k]     # k independent N(0,1) vectors
    return FrankWolfe.BlockVector(blocks)     # stack them as one block vector
end

# function used to generate a zero FrankWolfe.BlockVector, to test convex_feasibility_gradient variants
function zeros_blockvector(k::Integer, n::Integer; T=Float64)
    FrankWolfe.BlockVector([zeros(T, n) for _ in 1:k])
end