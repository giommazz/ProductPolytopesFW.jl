module BPCGProduct
using FrankWolfe

TODO: RIVEEDI LA FUNZIONE E COME ACCEDE AI DIVERSI BLOCCHI, POI FAI LO STESSO PER IL GRADIENTE

# Function to compute the pairwise distance objective
function objective_function(x, k, n)
    sum_dist = 0.0
    for i in 1:k    
        for j in i+1:k
            xi = x[(i-1)*n+1:i*n]
            xj = x[(j-1)*n+1:j*n]
            curr = sum((xi - xj).^2)
            println("°°°°°°°°°°°°°°°°°°°", x_i, x_j, curr)
            sum_dist += curr
        end
    end
    return 0.5 * sum_dist
end

# Gradient computation
function gradient!(storage, x, k, n)
    fill!(storage, 0.0)  # Reset storage
    for i in 1:k
        xi_index = (i-1)*n+1:i*n
        xi = x[xi_index]
        for j in 1:k
            if i != j
                xj_index = (j-1)*n+1:j*n
                xj = x[xj_index]
                storage[xi_index] .+= (xi - xj)
            end
        end
    end
end

# Setup Linear Minimization Oracles for the polytopes
function setup_lmos(n)
    # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
    lmo1 = FrankWolfe.ProbabilitySimplexOracle(1.0)
    lmo2 = FrankWolfe.ProbabilitySimplexOracle(12.0)
    lmo3 = FrankWolfe.ProbabilitySimplexOracle(5.0)
    lmo4 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))
    lmo5 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))
    lmo6 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), 3*ones(n))
    return (lmo1, lmo2, lmo3, lmo4, lmo5, lmo6)
end

# Initialize feasible points for each polytope
function initialize_points(lmos, n)
    starting_point =  tuple([FrankWolfe.compute_extreme_point(lmo, zeros(n)) for lmo in lmos]...)
    println("°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°", starting_point)
    return starting_point
end

# Function to run different optimization methods
function run_optimization_method(method, lmos, x0, k, n, target_tolerance, max_iterations)
    f = x -> objective_function(x, k, n)
    grad! = (storage, x) -> gradient!(storage, x, k, n)
    trajectories = []
    
    #     # Check (https://zib-iol.github.io/FrankWolfe.jl/dev/examples/docs_10_alternating_methods/): performing a full (simulatenous) BPCG update at each iteration, 
    #     #   by running `alternating_linear_minimization` with `blended_pairwise_conditional_gradient` inside
    x, v, primal, dual_gap, _, trajectory_data = FrankWolfe.alternating_linear_minimization(
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
    n = 2  # Dimension of each polytope
    k = 2   # Number of polytopes

    lmos = setup_lmos(n)
    # x0 = initialize_points(lmos, n)
    x0 = initialize_points(lmos, n)

    target_tolerance = 1e-6
    max_iterations = 10000

    # Run BPCG
    println("Running BPCG over the full polytope product...")
    bpcg_trajectories = run_optimization_method(FrankWolfe.blended_pairwise_conditional_gradient, lmos, x0, k, n, target_tolerance, max_iterations)

    # Run ALM with block-coordinate Frank-Wolfe
    println("Running ALM with block-coordinate Frank-Wolfe...")
    alm_trajectories = run_optimization_method(FrankWolfe.block_coordinate_frank_wolfe, lmos, x0, k, n, target_tolerance, max_iterations)

    println("Optimization completed.")
    #plot_trajectories([bpcg_trajectories, alm_trajectories], ["BPCG", "ALM"], xscalelog=true) # Uncomment and define plot_trajectories if visualization is needed
end

end # module BPCGProduct