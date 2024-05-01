module BPCGProduct
using FrankWolfe
using Plots
using Random

include(joinpath(dirname(pathof(FrankWolfe)), "../examples/plot_utils.jl"))

# Number of polytopes
k = 2
# Dimension of the subspace in which each polytope lies
n = 1000
# epsilon-optimality threshold
target_tolerance = 1e-9
# Number of FW iterations
max_iterations = 500000
# REMINDER FOR EXPERIMENTS: 
# 1) INCREASE precision
# 2) INCREASE N. OF MAX ITERATIONS SO THAT EVEN SUBLINEAR ALGO CONVERGES.
# 3) WITH THESE TWO, YOU CAN SHOW THAT TO REACH COVERGENCE BPCG IS FASTER 

function generate_rand_float_vector(lb=0, ub=100, seed=42)
    # Set the seed for reproducibility
    Random.seed!(seed)
    # Generate random Float64 in [a, b]
    return lb .+ (ub - lb) .* rand(Float64, n)
end

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

# Run BC (ALM setting) with vanilla FW
function run_FW(order, prod_lmo, x0, k, n, target_tolerance, max_iterations)   
    trajectories = []
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
# Run BC (ALM setting) with BPCG
function run_FW(order, update_step, prod_lmo, x0, k, n, target_tolerance, max_iterations)
    trajectories = []
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
# Run BPCG over full product LMO
function run_FW(prod_lmo, x0, k, n, target_tolerance, max_iterations)   
    trajectories = []
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

    # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
	# Setup Linear Minimization Oracles for the polytopes
    # lmo1 = FrankWolfe.ProbabilitySimplexOracle(11.0)                  # Probability simplex with given radius
    # bounds = generate_rand_float_vector()
    # lmo1 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))      # ℓ∞-norm ball 
    # lmo2 = FrankWolfe.ScaledBoundLInfNormBall(-bounds, bounds)        # scaled ℓ∞-norm ball
    # lmo1 = FrankWolfe.ScaledBoundL1NormBall(-ones(n), ones(n))        # ℓ₁-norm ball 
    # lmo2 = FrankWolfe.ScaledBoundL1NormBall(-bounds, bounds)          # scaled ℓ₁-norm ball 
    # lmo3 = FrankWolfe.ScaledBoundL1NormBall(-44*bounds, 5*bounds)     # scaled ℓ₁-norm ball 
    lmo1 = FrankWolfe.BirkhoffPolytopeLMO()
    lmo2 = FrankWolfe.BirkhoffPolytopeLMO()

    prod_lmo = create_product_lmo([lmo1, lmo2])
    println("°°°°°°°°°°°°°°°°", prod_lmo.lmos)
    x0 = find_starting_point(prod_lmo)

    trajectories = []

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

    #plot_trajectories([bc_fw_trajectories, bc_fw_trajectories, bpcg_trajectories], ["BC-FW", "BC-BPCG", "Full domain BPCG"], xscalelog=true)
end

end # module BPCGProduct