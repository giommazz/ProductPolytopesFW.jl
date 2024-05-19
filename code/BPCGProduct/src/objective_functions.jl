# `objective_functions.jl`

# Function to compute the pairwise distance objective: 1/2 ∑ᵢ₌₁ᵏ⁻¹∑ⱼ₌ᵢ₊₁ᵏ || xⁱ - xʲ ||₂²
function objective(config::Config, x::FrankWolfe.BlockVector)
    sum_dist = 0.0
    for i in 1:config.k  
        for j in i+1:config.k
            xi = x.blocks[i]
            xj = x.blocks[j]
            curr = sum((xi - xj).^2)
            sum_dist += curr
        end
    end
    return 0.5 * sum_dist
end

# Gradient computation for tuple of vectors
function gradient!(config::Config, storage::FrankWolfe.BlockVector, x::FrankWolfe.BlockVector)
    for i in 1:config.k
        sum_terms = zeros(config.n)
        for j in 1:config.k
            if i != j
                sum_terms .+= x.blocks[j]
            end
        end
        storage.blocks[i] .= 0.5 * ((config.k-1) * x.blocks[i] - sum_terms)
    end
end
