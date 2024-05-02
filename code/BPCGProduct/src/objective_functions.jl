# Function to compute the pairwise distance objective
function f(x)
    sum_dist = 0.0
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
function grad!(storage, x)
    for i = 1:k
        sum_terms = zeros(n)
        for j = 1:k
            if i != j
                sum_terms .+= x.blocks[j]
            end
        end
        storage.blocks[i] .= 0.5 * ((k-1) * x.blocks[i] - sum_terms)
    end
end