module BPCGProduct

using FrankWolfe
include("../examples/plot_utils.jl")

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
    lmo = FrankWolfe.ProductLMO(lmo1, lmo2)

    # Initialization
    p0 = [1; zeros(Int, n-1)]
    # Find initial feasible point for `lmo1`
    x01 = FrankWolfe.compute_extreme_point(lmo2, zeros(n)) #rand(n)
    # Find initial feasible point for `lmo2`
    x02 = FrankWolfe.compute_extreme_point(lmo1, zeros(n))
    
    # Run algorithm
    x, v, primal, dual_gap, _, trajectory_data = FrankWolfe.alternating_linear_minimization(    
        # TODO: check (https://zib-iol.github.io/FrankWolfe.jl/dev/examples/docs_10_alternating_methods/)
        #       I don't understand why there's a "FrankWolfe.block_coordinate_frank_wolfe" in the "Running Alternating Linear Minimization" example
        #       Also, on the same webgpage I read: "As an alternative to Block-Coordiante Frank-Wolfe (BCFW), one can also run alternating linear minimization with standard Frank-Wolfe algorithm. 
        #       These methods perform then the full (simulatenous) update at each iteration. In this example we also use FrankWolfe.away_frank_wolfe.". Should I use that?
        FrankWolfe.blended_pairwise_conditional_gradient,
        f,
        grad!,
        lmo, # TODO: should I use "lmos"? check (https://zib-iol.github.io/FrankWolfe.jl/dev/examples/docs_10_alternating_methods/)
        p0, 
        x0,
        update_order=FrankWolfe.FullUpdate(),
        epsilon=target_tolerance,
        max_iteration=k,
        #lazy=true, 
        line_search=FrankWolfe.Adaptive(L_est=L),  #FrankWolfe.Shortstep(one(T))
        print_iter=k / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    
    push!(trajectories, trajectory_data)
    plot_trajectories(trajectories, labels, xscalelog=true)
end





end # module