# `objective_functions.jl`

# Function to compute the pairwise distance objective
function convex_feasibility_objective(x::FrankWolfe.BlockVector)
    
    sum_dist = 0.0
    k = length(x.block_sizes)
    for i in 1:k  
        for j in i+1:k
            xi = x.blocks[i]
            xj = x.blocks[j]
            curr = sum((xi - xj).^2)
            sum_dist += curr
        end
    end
    # f(x) = 1/(2k) ∑ᵢ₌₁ᵏ⁻¹∑ⱼ₌ᵢ₊₁ᵏ ||xⁱ - xʲ||₂² = 1/(2k) [ (k-1) ∑ᵢ₌₁ᵏ||xⁱ||^2  -  2 ∑ᵢ₌₁ᵏ⁻¹∑ⱼ₌ᵢ₊₁ᵏ ⟨xⁱ, xʲ⟩ ]
    return 1.0 / (2k) * sum_dist
end

# Gradient computation for tuple of vectors
function convex_feasibility_gradient!(storage::FrankWolfe.BlockVector, x::FrankWolfe.BlockVector)
    
    k = length(x.block_sizes)
    n = Int(x.tot_size / k)
    for i in 1:k
        sum_terms = zeros(n)
        for j in 1:k
            if i != j
                sum_terms .+= x.blocks[j]
            end
        end
        # ∇ⁱf(x) = 1/k [ (k-1)xⁱ - ∑_{j ≠ i}xʲ ]
        storage.blocks[i] .= 1.0/k * ((k-1) * x.blocks[i] - sum_terms)
    end
end
