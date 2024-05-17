# `get_solutions.jl`
using BPCGProduct

config = Config("test/config.yml") # Use parameters from YAML file
filename = "intersecting_polytopes_n2_k10_v12-8-25-50-30-23-60-37-25-10_t20240516194119.jld2"

function get_solutions(config::Config, filename::String)

    # Load data and transform to Polyhedra.Polyhedron or JuMP.Model
    vertices, shifted_vertices = load_intersecting_polytopes(filename)

    # Initialize data structures
    lmo_list_jump_intersecting = Vector{FrankWolfe.MathOptLMO}()
    lmo_list_jump_nonintersecting = Vector{FrankWolfe.MathOptLMO}()
    
    lmo_list_cvxh_intersecting = Vector{FrankWolfe.ConvexHullOracle}()
    lmo_list_cvxh_nonintersecting = Vector{FrankWolfe.ConvexHullOracle}()    

    # Create JuMP.Model objects from each of the `vertices` entries
    for (v, vs) in zip(vertices, shifted_vertices)    
        # Create objects of type Polyhedra.Polyhedron from `v` and `vs`
        poly_intersecting = polytope(vs)
        poly_nonintersecting = polytope(v)
        # Create objects of type JuMP.Model from `poly` and `poly_shifted`
        jump_poly_intersecting = polyhedra_to_jump(config, poly_intersecting)
        jump_poly_nonintersecting = polyhedra_to_jump(config, poly_nonintersecting)
        lmo_intersecting = FrankWolfe.MathOptLMO(jump_poly_intersecting.moi_backend)
        lmo_nonintersecting = FrankWolfe.MathOptLMO(jump_poly_nonintersecting.moi_backend)
        push!(lmo_list_jump_intersecting, lmo_intersecting)
        push!(lmo_list_jump_nonintersecting, lmo_nonintersecting)
    end

    
    
        
        
    end
    for v in vertices
        lmo = FrankWolfe.ConvexHullOracle(v)
        push!(lmo_list, lmo)
    end

    # TODO: CREATE LMO FROM 
    #   -   JUMP MODEL: something like lmo = FrankWolfe.MathOptLMO(model.moi_backend) OR 
    #   -   VERTICES: something like lmo = FrankWolfe.ConvexHullOracle([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    lmo_list = []#TODO: put all lmos in a list
    # Find possible subsets of size `config.k`
    lmo_products = unique_combinations(lmo_list, config)
    
    for lmos in lmo_products
        prod_lmo = create_product_lmo(config, lmos)
        println("\n\n\n---------------------------------------------------------")
        println("---------------------------------------------------------")
        println("LMOs: ", [typeof(prod_lmo.lmos[i]) for i in 1:config.k])
        println("---------------------------------------------------------")
        println("---------------------------------------------------------")
        x0 = find_starting_point(config, prod_lmo)
        # Block-coordinate BPCG with CyclicUpdate
        println("\n\n\n ----------> Full Block-coordinate BPCG")
        trajectories = run_FW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo, x0)
    end
end
# 4. have a function that does 1, 2, 3 and then executes fw run_fw_variants
#         run_fw_variants(lmos, fw_variant) # takes lmo as inputs, and a fw variant as dict
#         run_fw_variants(lmos, fw_varians) # takes lmo as inputs, and multiple fw variants as dict
#         run_fw_variants(vertices, fw_variant) # takes polytope vertices as inputs, and a fw variant as dict
#         run_fw_variants(vertices, fw_variants) # takes polytope vertices as inputs, and multiple fw variants as dict


# Call the main function to run your experiments
# function main(config::Config)

#     println("k: $(config.k), n: $(config.n), target_tolerance: $(config.target_tolerance), max_iterations: $(config.max_iterations), seed: $(config.seed)")

#     # Polytopes LMOs (https://github.com/ZIB-IOL/FrankWolfe.jl/blob/97a599c029a054aab6a5574d9bed8d48e0f9fb01/src/polytope_oracles.jl#L4)
#     # Setup Linear Minimization Oracles for the polytopes
#     lmo_probsmplx_1 = FrankWolfe.ProbabilitySimplexOracle(5.0)                    # Probability simplex with given radius
#     lmo_probsmplx_2 = FrankWolfe.ProbabilitySimplexOracle(1.0)                      # Probability simplex with given radius
#     lmo_probsmplx_3 = FrankWolfe.ProbabilitySimplexOracle(50.0)                     # Probability simplex with given radius
#     bounds = generate_rand_float_vector(config)
#     lmo_infnormball_1 = FrankWolfe.ScaledBoundLInfNormBall(-ones(config.n), ones(config.n))     # ℓ∞-norm ball
#     lmo_infnormball_2 = FrankWolfe.ScaledBoundLInfNormBall(-bounds, bounds)                     # scaled ℓ∞-norm ball
#     lmo_infnormball_3 = FrankWolfe.ScaledBoundLInfNormBall(bounds, 5*bounds)                     # scaled ℓ∞-norm ball
#     lmo_onenormball_1 = FrankWolfe.ScaledBoundL1NormBall(6*ones(config.n), 15*ones(config.n))       # ℓ₁-norm ball
#     lmo_onenormball_2 = FrankWolfe.ScaledBoundL1NormBall(bounds, 50*bounds)                      # scaled ℓ₁-norm ball
#     lmo_onenormball_3 = FrankWolfe.ScaledBoundL1NormBall(-bounds, bounds)                       # scaled ℓ₁-norm ball 
# end

# main(config)
