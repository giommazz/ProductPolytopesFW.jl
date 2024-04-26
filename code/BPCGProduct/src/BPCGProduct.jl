module BPCGProduct
using FrankWolfe


function grad!(storage, x)
    @. storage = zero(x)
end

function run_bpcg_fulldomain()
    
    println("Running BPCG over the full polytope product...")
    
    # Parameters
    n = 20
    k = 10000
    f(x) = 0
    target_tolerance = 1e-9
    trajectories = [];

    # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
    lmo1 = FrankWolfe.ProbabilitySimplexOracle(1.0)
    lmo2 = FrankWolfe.ProbabilitySimplexOracle(12.0)
    lmo3 = FrankWolfe.ProbabilitySimplexOracle(5.0)
    lmo4 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))
    lmo5 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))
    lmo6 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), 3*ones(n))
    lmos = (lmo1, lmo2, lmo3, lmo4, lmo5, lmo6)
    
    # Find initial feasible points
    x01 = FrankWolfe.compute_extreme_point(lmo1, zeros(n))
    x02 = FrankWolfe.compute_extreme_point(lmo2, zeros(n))
    x03 = FrankWolfe.compute_extreme_point(lmo3, zeros(n))
    x04 = FrankWolfe.compute_extreme_point(lmo4, zeros(n))
    x05 = FrankWolfe.compute_extreme_point(lmo5, zeros(n))
    x06 = FrankWolfe.compute_extreme_point(lmo6, zeros(n))
    x0 = (x01, x02, x03, x04, x05, x06)

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


function run_alm()
    
    println("Running BPCG over the full polytope product...")
    
    # Parameters
    n = 20
    k = 10000
    f(x) = 0
    target_tolerance = 1e-9
    trajectories = [];

    # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
    lmo1 = FrankWolfe.ProbabilitySimplexOracle(1.0)
    lmo2 = FrankWolfe.ProbabilitySimplexOracle(12.0)
    lmo3 = FrankWolfe.ProbabilitySimplexOracle(5.0)
    lmo4 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))
    lmo5 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), zeros(n))
    lmo6 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), 3*ones(n))
    lmos = (lmo1, lmo2, lmo3, lmo4, lmo5, lmo6)
    
    # Find initial feasible points
    x01 = FrankWolfe.compute_extreme_point(lmo1, zeros(n))
    x02 = FrankWolfe.compute_extreme_point(lmo2, zeros(n))
    x03 = FrankWolfe.compute_extreme_point(lmo3, zeros(n))
    x04 = FrankWolfe.compute_extreme_point(lmo4, zeros(n))
    x05 = FrankWolfe.compute_extreme_point(lmo5, zeros(n))
    x06 = FrankWolfe.compute_extreme_point(lmo6, zeros(n))
    x0 = (x01, x02, x03, x04, x05, x06)

    x, v, primal, dual_gap, _, trajectory_data = FrankWolfe.alternating_linear_minimization(
        FrankWolfe.block_coordinate_frank_wolfe,
        f,
        grad!,
        lmos,
        x0,
        epsilon=target_tolerance,
        max_iteration=k,
        update_order=FrankWolfe.CyclicUpdate(),
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=k / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    )
    push!(trajectories, trajectory_data)
end



end # module