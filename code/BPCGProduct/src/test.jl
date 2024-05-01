module test
using FrankWolfe

k = 2
n = 5000
target_tolerance = 1e-9
max_iterations = 500000
# REMINDER FOR EXPERIMENTS: 
# 1) INCREASE precision
# 2) INCREASE N. OF MAX ITERATIONS SO THAT EVEN SUBLINEAR ALGO CONVERGES.
# 3) WITH THESE TWO, YOU CAN SHOW THAT TO REACH COVERGENCE BPCG IS FASTER 


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

function create_product_lmo(lmos_list)
    # Check if length of LMO list matches `k`
    if length(lmos_list) != k
        error("The number of LMOs provided ($(length(lmos_list))) does not match the expected number ($k).")
    end
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmos_list)
    # Create and return a ProductLMO object
    return FrankWolfe.ProductLMO(lmos_tuple)
end

# Find starting point `x0` over the product of different LMOs
function find_starting_point(prod_lmo)    
    # Prepare datafor `x0`, which is of type `FrankWolfe.BlockVector{Float64, Vector{Float64}, Tuple{Int64}}`
    # 1) Compute extreme points for each LMO in the product
    extreme_points = [FrankWolfe.compute_extreme_point(lmo, zeros(n)) for lmo in prod_lmo.lmos]    
    # 2) Generate a vector of length `k`, where each entry is a tuple `(n,)`
    block_sizes = fill((n,), length(prod_lmo.lmos))
    # 3) Compute total size of `x0`
    total_size = sum([size[1] for size in block_sizes]) 

    # Instantiate `x0`
    return FrankWolfe.BlockVector(extreme_points, block_sizes, total_size)
end

# Run ALM with vanilla FW
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
# Run ALM with BPCG
function run_FW(order, update_step, prod_lmo, x0, k, n, target_tolerance, max_iterations)
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
        update_step=update_step,
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

function main()

	# Setup Linear Minimization Oracles for the polytopes
    lmo1 = FrankWolfe.ProbabilitySimplexOracle(1.0)
    lmo2 = FrankWolfe.ProbabilitySimplexOracle(123.0)
    # lmo3 = FrankWolfe.ProbabilitySimplexOracle(123.0)
    # lmo4 = FrankWolfe.ProbabilitySimplexOracle(123.0)
    # lmo5 = FrankWolfe.ProbabilitySimplexOracle(123.0)
    # lmo6 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))
    # lmo7 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))
    # lmo8 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))
    # lmo9 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))
    # lmo10 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))

    prod_lmo = create_product_lmo([lmo1, lmo2])
    x0 = find_starting_point(prod_lmo)

    # Block-coordinate vanilla FW
    println("\n\n\n---------------------------------------------------------")
    println("Run Alternating Linear Minimization (block-coordinate style) with vanilla FW")
    println("---------------------------------------------------------")
    bc_fw_trajectories = run_FW(FrankWolfe.CyclicUpdate(), prod_lmo, x0, k, n, target_tolerance, max_iterations)
    
    # Block-coordinate BPCG
    println("\n\n\n---------------------------------------------------------")
    println("Run Alternating Linear Minimization (block-coordinate style) with BPCG)")
    println("---------------------------------------------------------")
    bc_bpcg_trajectories = run_FW(FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo, x0, k, n, target_tolerance, max_iterations)
    
    # BPCG over full product LMO
    println("\n\n\n---------------------------------------------------------")
    println("Run Blended Pairwise Conditional Gradient over the full product LMO")
    println("---------------------------------------------------------")
    bpcg_trajectories = run_FW(prod_lmo, x0, k, n, target_tolerance, max_iterations)
end

end # module BPCGProduct