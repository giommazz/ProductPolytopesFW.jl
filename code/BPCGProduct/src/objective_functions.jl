# `objective_functions.jl`

"""
Function to compute form 1 of the objective function
f(x)  =     1/(2k) ∑_{1 ≤ i < j ≤ k} ||xⁱ - xʲ||₂²           [form 1]
      =     1/(2k) [ (k-1) ∑ᵢ₌₁ᵏ||xⁱ||²  -  2 ∑_{1 ≤ i < j ≤ k} ⟨xⁱ, xʲ⟩ ]       [form 2.a]
      =     1/(2k) k∑ᵢ₌₁ᵏ||xⁱ||² - ||∑_{1 ≤ i < j ≤ k} xⁱ||²                    [form 2.b]
      =     1/(2k)<x, Mₖ x>, where Mₖ = kI - 𝟏𝟏ᵀ, i.e., a (kn X kn) matrix with (k-1) on the diagonal and -1 elsewhere           [form 3]


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
function convex_feasibility_objective_v2b(x::FrankWolfe.BlockVector)
    
    # number of blocks
    k = length(x.block_sizes)
    # dimension of each block
    n = length(x.blocks[1])
    # accumulator for Σ‖xᵢ‖²
    s_norm = 0.0
    # Σ xᵢ
    s_vec = Vector{eltype(x)}(undef, n)
    # zero‑initialise once
    fill!(s_vec, 0)

    # fast, bounds‑check‑free loop
    @inbounds @simd for blk in x.blocks
        # BLAS level‑1
        s_norm += dot(blk, blk)
        # s_vec += blk   (BLAS)
        LinearAlgebra.BLAS.axpy!(1.0, blk, s_vec)
    end
    return 0.5/k * (k*s_norm - dot(s_vec, s_vec))
end


# Gradient computation for tuple of vectors
# ∇ⁱf(x) = 1/k [ (k-1)xⁱ - ∑_{j ≠ i}xʲ ]
function convex_feasibility_gradient!(storage::FrankWolfe.BlockVector, x::FrankWolfe.BlockVector)
    
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

# TODO: THIS FUNCTION IS JUST USED TO TEST THE OBJECTIVE FUNCTIONS. DELETE WHEN DONE
function random_blockvector(k::Integer, n::Integer; rng = Random.default_rng())
    blocks = [randn(rng, n) for _ in 1:k]     # k independent N(0,1) vectors
    return FrankWolfe.BlockVector(blocks)     # stack them as one block vector
end
