module BPCGProduct
using FrankWolfe

k = 2
n = 2


# Function to compute the pairwise distance objective
function f(x)
    sum_dist = 0.0
    for i in 1:k    
        for j in i+1:k
            xi = x[i]
            xj = x[j]
            curr = sum((xi - xj).^2)
            sum_dist += curr
        end
    end
    return 0.5 * sum_dist
end

# Gradient computation for tuple of vectors
function grad!(storage, x)
    # Initialize or reset 'storage' to zero vectors
    storage = [zeros(Float64, n) for _ in 1:k]
    println("°°°°°°°° storage: ", storage, typeof(storage))
    
    println("°°°°°°°°x: ", x, typeof(x))
    for i = 1:k
        println("°°°°°°°° xᵢ: ", x[i], typeof(x))
        sum_terms = zeros(n)
        for j = 1:k
            if i != j
                println("°°°°°°°° xⱼ", x[j], typeof(x))
                sum_terms .+= x[j]
            end
        end
        storage[i] .= 0.5 * ((k-1) * x[i] - sum_terms)
    end
end

function grad!(storage, x)
    # Initialize or reset 'storage' to zero vectors
    storage = [zeros(Float64, n) for _ in 1:k]
    println("°°°°°°°° storage: ", storage, typeof(storage))
    
    println("°°°°°°°°x: ", x, typeof(x))
    for i = 1:k
        println("°°°°°°°° xᵢ: ", x[i], typeof(x))
        sum_terms = zeros(n)
        for j = 1:k
            if i != j
                println("°°°°°°°° xⱼ", x[j], typeof(x))
                sum_terms .+= x[j]
            end
        end
        storage[i] .= 0.5 * ((k-1) * x[i] - sum_terms)
    end
end




# Setup Linear Minimization Oracles for the polytopes
function setup_lmos(n)
    # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
    lmo1 = FrankWolfe.ProbabilitySimplexOracle(1.0)
    lmo2 = FrankWolfe.ProbabilitySimplexOracle(12.0)
    #lmo3 = FrankWolfe.ProbabilitySimplexOracle(5.0)
    #lmo4 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))
    #lmo5 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))
    #lmo6 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), 3*ones(n))
    return (lmo1, lmo2#, lmo3, lmo4, lmo5, lmo6
    )
end

# Initialize feasible points for each polytope
function initialize_points(lmos, n)
    return  tuple([FrankWolfe.compute_extreme_point(lmo, zeros(n)) for lmo in lmos]...)
end

# Function to run different optimization methods
function run_optimization_method(method, lmos, x0, k, n, target_tolerance, max_iterations)
    
    
    #storage = [zeros(Float64, n) for _ in 1:k]
    #f = x -> objective_function(x, k, n)
    #grad! = (storage, x) -> gradient!(storage, x, k, n)
    
    trajectories = []
    println("°°°°°°°° x₀: ", x0, typeof(x0))

    #     # Check (https://zib-iol.github.io/FrankWolfe.jl/dev/examples/docs_10_alternating_methods/): performing a full (simulatenous) BPCG update at each iteration, 
    #     #   by running `alternating_linear_minimization` with `blended_pairwise_conditional_gradient` inside
    xx, v, primal, dual_gap, _, trajectory_data = FrankWolfe.alternating_linear_minimization(
        method,
        f,
        grad!,
        lmos,
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
    #plot_trajectories(trajectories, labels, xscalelog=true)
    
    return trajectories
end

# Main execution function
function main()

    lmos = setup_lmos(n)
    x0 = initialize_points(lmos, n)

    target_tolerance = 1e-6
    max_iterations = 10000

    # Run BPCG
    println("Running BPCG over the full polytope product...")
    bpcg_trajectories = run_optimization_method(FrankWolfe.blended_pairwise_conditional_gradient, lmos, x0, k, n, target_tolerance, max_iterations)

    # Run ALM with block-coordinate Frank-Wolfe
    #println("Running ALM with block-coordinate Frank-Wolfe...")
    #alm_trajectories = run_optimization_method(FrankWolfe.block_coordinate_frank_wolfe, lmos, x0, k, n, target_tolerance, max_iterations)

    println("Optimization completed.")
    #plot_trajectories([bpcg_trajectories, alm_trajectories], ["BPCG", "ALM"], xscalelog=true) # Uncomment and define plot_trajectories if visualization is needed
end

end # module BPCGProduct