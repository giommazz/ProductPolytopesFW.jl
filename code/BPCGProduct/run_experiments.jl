# Assuming your module's code is in 'src/BPCGProduct.jl'
include("src/BPCGProduct.jl")  # Include the module
using .BPCGProduct  # Use the module

# Call the main function to run your experiments
function main()
        
    println("k = $k, n = $n, seed = $seed")
    
    lmos = [1, 2, 3]
    println(unique_combinations(lmos))
    readline()


    # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
    # Setup Linear Minimization Oracles for the polytopes
    lmo_probsmplx_1 = FrankWolfe.ProbabilitySimplexOracle(200.0)                    # Probability simplex with given radius
    lmo_probsmplx_2 = FrankWolfe.ProbabilitySimplexOracle(1.0)                      # Probability simplex with given radius
    lmo_probsmplx_3 = FrankWolfe.ProbabilitySimplexOracle(50.0)                     # Probability simplex with given radius
    bounds = generate_rand_float_vector()
    lmo_infnormball_1 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))       # ℓ∞-norm ball
    lmo_infnormball_2 = FrankWolfe.ScaledBoundLInfNormBall(-bounds, bounds)         # scaled ℓ∞-norm ball
    lmo_onenormball_1 = FrankWolfe.ScaledBoundL1NormBall(-ones(n), ones(n))         # ℓ₁-norm ball
    lmo_onenormball_2 = FrankWolfe.ScaledBoundL1NormBall(bounds, 4*bounds)          # scaled ℓ₁-norm ball
    lmo_onenormball_3 = FrankWolfe.ScaledBoundL1NormBall(-bounds, bounds)           # scaled ℓ₁-norm ball 
    
    lmos = [lmo_probsmplx_1, lmo_probsmplx_2, lmo_probsmplx_3, lmo_infnormball_1, lmo_infnormball_2, lmo_onenormball_1, lmo_onenormball_2, lmo_onenormball_3]
    
    
    prod_lmo = create_product_lmo([lmo1, lmo2])
    println("LMOs: ", prod_lmo.lmos)
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

    #plot_trajectories([bc_fw_trajectories, bc_fw_trajectories, bpcg_trajectories], ["BC-FW", "BC-BPCG", "Full domain BPCG"], xscalelog=true)
end

