# `run_experiments.jl`
using BPCGProduct
using BPCGProduct: FrankWolfe

config = Config("test/config.yml") # Use parameters from YAML file

# Call the main function to run your experiments
function main(config)

    println("k: $(config.k), n: $(config.n), target_tolerance: $(config.target_tolerance), max_iterations: $(config.max_iterations), seed: $(config.seed)")

    # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
    # Setup Linear Minimization Oracles for the polytopes
    lmo_probsmplx_1 = FrankWolfe.ProbabilitySimplexOracle(200.0)                    # Probability simplex with given radius
    lmo_probsmplx_2 = FrankWolfe.ProbabilitySimplexOracle(1.0)                      # Probability simplex with given radius
    lmo_probsmplx_3 = FrankWolfe.ProbabilitySimplexOracle(50.0)                     # Probability simplex with given radius
    bounds = generate_rand_float_vector(config)
    lmo_infnormball_1 = FrankWolfe.ScaledBoundLInfNormBall(-ones(config.n), ones(config.n))     # ℓ∞-norm ball
    lmo_infnormball_2 = FrankWolfe.ScaledBoundLInfNormBall(-bounds, bounds)                     # scaled ℓ∞-norm ball
    lmo_infnormball_3 = FrankWolfe.ScaledBoundLInfNormBall(bounds, 5*bounds)                     # scaled ℓ∞-norm ball
    lmo_onenormball_1 = FrankWolfe.ScaledBoundL1NormBall(-ones(config.n), ones(config.n))       # ℓ₁-norm ball
    lmo_onenormball_2 = FrankWolfe.ScaledBoundL1NormBall(bounds, 50*bounds)                      # scaled ℓ₁-norm ball
    lmo_onenormball_3 = FrankWolfe.ScaledBoundL1NormBall(-bounds, bounds)                       # scaled ℓ₁-norm ball 
    
    #lmos = [lmo_probsmplx_1, lmo_probsmplx_2, lmo_probsmplx_3, lmo_infnormball_1, lmo_infnormball_2, lmo_onenormball_1, lmo_onenormball_2, lmo_onenormball_3]
    lmo_list = [lmo_probsmplx_1, lmo_probsmplx_3, lmo_infnormball_3, lmo_onenormball_2]
    lmo_products = unique_combinations(lmo_list, config)
    
    for lmos in lmo_products
        
        prod_lmo = create_product_lmo(lmos, config)
        println("\n\n\n---------------------------------------------------------")
        println("---------------------------------------------------------")
        println("LMOs: ", [typeof(prod_lmo.lmos[i]) for i in 1:config.k])
        println("---------------------------------------------------------")
        println("---------------------------------------------------------")
        x0 = find_starting_point(prod_lmo, config)

        # Block-coordinate vanilla FW
        println("\n\n\n ----------> Block-coordinate vanilla FW")
        bc_fw_trajectories = run_FW(FrankWolfe.CyclicUpdate(), prod_lmo, x0, config)

        # Block-coordinate BPCG
        println("\n\n\n ----------> Block-coordinate BPCG")
        bc_bpcg_trajectories = run_FW(FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo, x0, config)
        
        # TODO: TRY THIS! THIS SHOULD BE THE SAME AS CALLING FW.BPCG. IF NOT, SOMETHING IS WRONG
        println("\n\n\n ----------> Block-coordinate BPCG over FULL product (i.e., BC with full BPCG updates)")
        bpcg_trajectories = run_FW(FrankWolfe.FullUpdate(), FrankWolfe.BPCGStep(), prod_lmo, x0, config)

        # BPCG over full product LMO
        println("\n\n\n ----------> BPCG over full product LMO")
        bpcg_trajectories = run_FW(prod_lmo, x0, config)

        #plot_trajectories([bc_fw_trajectories, bc_fw_trajectories, bpcg_trajectories], ["BC-FW", "BC-BPCG", "Full domain BPCG"], xscalelog=true)
    end
end

main(config)
