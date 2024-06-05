# `objective_functions.jl`

# Function to compute the pairwise distance objective: 1/2 ∑ᵢ₌₁ᵏ⁻¹∑ⱼ₌ᵢ₊₁ᵏ || xⁱ - xʲ ||₂²
function objective(x::FrankWolfe.BlockVector)
    
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
    return 0.5 * sum_dist
end
# (Multiple dispatch)
function objective(x::FrankWolfe.BlockVector)
    
    TODO: FAI FUNZIONE PER LMO SEMPLICE CON NORMA

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
    return 0.5 * sum_dist
end

# Gradient computation for tuple of vectors
function gradient!(storage::FrankWolfe.BlockVector, x::FrankWolfe.BlockVector)
    
    k = length(x.block_sizes)
    n = Int(x.tot_size / k)
    for i in 1:k
        sum_terms = zeros(n)
        for j in 1:k
            if i != j
                sum_terms .+= x.blocks[j]
            end
        end
        storage.blocks[i] .= 0.5 * ((k-1) * x.blocks[i] - sum_terms)
    end
end
