# `get_solutions.jl`
using BPCGProduct

config = Config("test/config.yml") # Use parameters from YAML file
filename = "intersecting_polytopes_n2_k10_v12-8-25-50-30-23-60-37-25-10_t20240516194119.jld2"

# Load data and transform to Polyhedra.Polyhedron or JuMP.Model
vertices, shifted_vertices = load_intersecting_polytopes(filename)
polytopes_polyhedra = Vector{Polyhedron}()
polytopes_jump = Vector{Model}()
for v in vertices
    poly = polytope(v)
    push!(polytopes_polyhedra, poly)
    push!(polytopes_polyhedra, polyhedra_to_jump(config, poly))
end

# TODO: CREATE FRANKWOLFE LMOS FROM polytopes_polyhedra and from polytopes_jump


# 1. get instances
# 2. define LMOs
# 3. lmo list and lmo products
# 4. have a function that does 1, 2, 3 and then executes fw run_fw_variants
#         run_fw_variants(lmos, fw_variant) # takes lmo as inputs, and a fw variant as dict
#         run_fw_variants(lmos, fw_varians) # takes lmo as inputs, and multiple fw variants as dict
#         run_fw_variants(vertices, fw_variant) # takes polytope vertices as inputs, and a fw variant as dict
#         run_fw_variants(vertices, fw_variants) # takes polytope vertices as inputs, and multiple fw variants as dict


# Call the main function to run your experiments
function main(config::Config)

    println("k: $(config.k), n: $(config.n), target_tolerance: $(config.target_tolerance), max_iterations: $(config.max_iterations), seed: $(config.seed)")

    # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
    # Setup Linear Minimization Oracles for the polytopes
    lmo_probsmplx_1 = FrankWolfe.ProbabilitySimplexOracle(5.0)                    # Probability simplex with given radius
    lmo_probsmplx_2 = FrankWolfe.ProbabilitySimplexOracle(1.0)                      # Probability simplex with given radius
    lmo_probsmplx_3 = FrankWolfe.ProbabilitySimplexOracle(50.0)                     # Probability simplex with given radius
    bounds = generate_rand_float_vector(config)
    lmo_infnormball_1 = FrankWolfe.ScaledBoundLInfNormBall(-ones(config.n), ones(config.n))     # ℓ∞-norm ball
    lmo_infnormball_2 = FrankWolfe.ScaledBoundLInfNormBall(-bounds, bounds)                     # scaled ℓ∞-norm ball
    lmo_infnormball_3 = FrankWolfe.ScaledBoundLInfNormBall(bounds, 5*bounds)                     # scaled ℓ∞-norm ball
    lmo_onenormball_1 = FrankWolfe.ScaledBoundL1NormBall(6*ones(config.n), 15*ones(config.n))       # ℓ₁-norm ball
    lmo_onenormball_2 = FrankWolfe.ScaledBoundL1NormBall(bounds, 50*bounds)                      # scaled ℓ₁-norm ball
    lmo_onenormball_3 = FrankWolfe.ScaledBoundL1NormBall(-bounds, bounds)                       # scaled ℓ₁-norm ball 

    #lmo_list = [lmo_probsmplx_1, lmo_probsmplx_3, lmo_infnormball_1, lmo_infnormball_3, lmo_onenormball_2, lmo_onenormball_3]
    lmo_list = [lmo_probsmplx_1, lmo_onenormball_1]
    lmo_products = unique_combinations(lmo_list, config)
    
    for lmos in lmo_products
        
        prod_lmo = create_product_lmo(config, lmos)
        println("\n\n\n---------------------------------------------------------")
        println("---------------------------------------------------------")
        println("LMOs: ", [typeof(prod_lmo.lmos[i]) for i in 1:config.k])
        println("---------------------------------------------------------")
        println("---------------------------------------------------------")
        x0 = find_starting_point(config, prod_lmo)

        # Block-coordinate vanilla FW
        println("\n\n\n ----------> Cyclic Block-coordinate vanilla FW")
        bc_fw_trajectories = run_FW(config, FrankWolfe.CyclicUpdate(), prod_lmo, x0)

        # Block-coordinate BPCG
        println("\n\n\n ----------> Cyclic Block-coordinate BPCG")
        bc_bpcg_cyclic_trajectories = run_FW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo, x0)
        
        # TODO: TRY THIS! THIS SHOULD BE THE SAME AS CALLING FW.BPCG. IF NOT, SOMETHING IS WRONG
        println("\n\n\n ----------> Full Block-coordinate BPCG")
        bc_bpcg_full_trajectories = run_FW(config, FrankWolfe.FullUpdate(), FrankWolfe.BPCGStep(), prod_lmo, x0)

        # BPCG over full product LMO
        println("\n\n\n ----------> BPCG")
        bpcg_trajectories = run_FW(config, prod_lmo, x0)

        println("\n\n\n ----------> AP")
        altproj_trajectories = run_FW(config, prod_lmo, x0, true)

        #plot_trajectories([bc_fw_trajectories, bc_fw_trajectories, bpcg_trajectories], ["BC-FW", "BC-BPCG", "Full domain BPCG"], xscalelog=true)
    end
end

main(config)
