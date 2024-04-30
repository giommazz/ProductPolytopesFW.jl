module test
using FrankWolfe

k = 3
n = 2
target_tolerance = 1e-6
max_iterations = 10000


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


# TODO: SEE BLOCKVECTOR STRUCTURE AND PROD_LMO STRUCTURE

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

# Run ALM
function run_FW(order, prod_lmo, x0, k, n, target_tolerance, max_iterations)   
    trajectories = []

    # Check (https://zib-iol.github.io/FrankWolfe.jl/dev/examples/docs_10_alternating_methods/): performing a full (simulatenous) BPCG update at each iteration, 
    # 	by running `alternating_linear_minimization` with `blended_pairwise_conditional_gradient` inside
    x, v, primal, dual_gap, trajectory_data = FrankWolfe.block_coordinate_frank_wolfe(
        f,
        grad!,
        prod_lmo,
        x0,
        update_order=order,
        epsilon=target_tolerance,
        max_iteration=max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
# Run BPCG
function run_FW(prod_lmo, x0, k, n, target_tolerance, max_iterations)   
    trajectories = []
    # Check (https://zib-iol.github.io/FrankWolfe.jl/dev/examples/docs_10_alternating_methods/): performing a full (simulatenous) BPCG update at each iteration, 
    # 	by running `alternating_linear_minimization` with `blended_pairwise_conditional_gradient` inside
    x, v, primal, dual_gap, _, trajectory_data = FrankWolfe.blended_pairwise_conditional_gradient(
        f,
        grad!,
        prod_lmo,
        x0,
        epsilon=target_tolerance,
        max_iteration=max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
	
# Main execution function
function main()

	# Setup Linear Minimization Oracles for the polytopes
    lmo1 = FrankWolfe.ProbabilitySimplexOracle(1.0)
    lmo2 = FrankWolfe.ProbabilitySimplexOracle(12.0)
    lmo3 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))

    lmos = (lmo1, lmo2, lmo3)
    prod_lmo = FrankWolfe.ProductLMO(lmos)

    x0 = FrankWolfe.BlockVector([-ones(n), [i == 1 ? 1 : 0 for i in 1:n]], fill((n,), k), k * n)
    # src/alternating_linear_minimization FrankWolfe.jl
    #x0 = compute_extreme_point(FrankWolfe.ProductLMO(lmos), tuple(fill(zeros(n), k)...))
    #x0 = tuple([FrankWolfe.compute_extreme_point(lmo, zeros(n)) for lmo in lmos]...)

    # Run BPCG: FrankWolfe.blended_pairwise_conditional_gradient, 
    bpcg_trajectories = run_FW(prod_lmo, x0, k, n, target_tolerance, max_iterations)

    # Run ALM with block-coordinate Frank-Wolfe: FrankWolfe.block_coordinate_frank_wolfe
    alm_trajectories = run_FW(FrankWolfe.CyclicUpdate(), prod_lmo, x0, k, n, target_tolerance, max_iterations)

end

end # module BPCGProduct