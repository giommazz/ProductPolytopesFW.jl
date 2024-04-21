module BPCGProduct
using FrankWolfe


function grad!(storage, x)
    @. storage = zero(x)
end

# Example function that utilizes FrankWolfe
function runExperiment()
    
    println("Running experiment with Frank Wolfe...")
    
    # Parameters
    n = 20
    k = 10000
    f(x) = 0
    target_tolerance = 1e-6
    trajectories = [];

    # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
    lmo1 = FrankWolfe.ProbabilitySimplexOracle(1.0)
    lmo2 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))
    lmos = (lmo1, lmo2)

    # Initialization
    # Find initial feasible point for `lmo1`
    x01 = FrankWolfe.compute_extreme_point(lmo2, zeros(n)) #rand(n)
    # Find initial feasible point for `lmo2`
    x02 = FrankWolfe.compute_extreme_point(lmo1, zeros(n))
    x0 = (x01, x02)

    # Run algorithm
    x, v, primal, dual_gap, _, trajectory_data = FrankWolfe.alternating_linear_minimization(    
        # Check (https://zib-iol.github.io/FrankWolfe.jl/dev/examples/docs_10_alternating_methods/): performing a full (simulatenous) BPCG update at each iteration, 
        #   by running `alternating_linear_minimization` with `blended_pairwise_conditional_gradient` inside
        FrankWolfe.blended_pairwise_conditional_gradient,
        f,
        grad!,
        lmos,
        x0,  
        epsilon=target_tolerance,
        max_iteration=k,
        #lazy=true, 
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=k / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    
    push!(trajectories, trajectory_data)
    #plot_trajectories(trajectories, labels, xscalelog=true)
end

end # module