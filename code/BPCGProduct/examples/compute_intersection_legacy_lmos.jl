# `compute_intersection_legacy_lmos.jl`
using ProductPolytopesAFW
using FrankWolfe

config = Config("examples/config.yml") # Use parameters from YAML file

function lmo_name(lmo_list::Vector{FrankWolfe.LinearMinimizationOracle})
    if lmo_list == [lmo_probsmplx_1, lmo_unitsmplx] return "prob_unit_i", false
    elseif lmo_list == [lmo_probsmplx_2, lmo_unitsmplx] return "prob_unit_ni", true
    
    elseif lmo_list == [lmo_unitsmplx, lmo_onenormball_1] return "intersecting_l1_unit_i", false
    elseif lmo_list == [lmo_unitsmplx, lmo_onenormball_2] return "intersecting_l1_unit_i", false
    elseif lmo_list == [lmo_unitsmplx, lmo_onenormball_3] return "intersecting_l1_unit_ni", true
    
    elseif lmo_list == [lmo_unitsmplx, lmo_infnormball_1] return "intersecting_linf_unit_i", false
    elseif lmo_list == [lmo_unitsmplx, lmo_infnormball_2] return "intersecting_linf_unit_i", false
    elseif lmo_list == [lmo_unitsmplx, lmo_infnormball_4] return "intersecting_linf_unit_i", false
    elseif lmo_list == [lmo_unitsmplx, lmo_infnormball_3] return "intersecting_linf_unit_ni", true
    
    elseif lmo_list == [lmo_infnormball_2, lmo_onenormball_4] return "intersecting_linf_l1_i", false
    elseif lmo_list == [lmo_infnormball_2, lmo_onenormball_3] return "intersecting_linf_l1_i", false
    elseif lmo_list == [lmo_infnormball_1, lmo_onenormball_2] return "intersecting_linf_l1_ni", true
    end
end

# Call the main function to run your experiments
function main(config::Config, lmo_list::Vector{FrankWolfe.LinearMinimizationOracle}, ni_flag::Bool, basename::String)

    trajectories = []
        
    primal, _ = compute_distance(config, lmo_list)

    prod_lmo = create_product_lmo(config, lmo_list)
    println("\n\n\n---------------------------------------------------------")
    println("LMOs in ProductLMO: ")
    for i in [typeof(prod_lmo.lmos[i]) for i in 1:config.k] println("\t$i") end
    println("---------------------------------------------------------")

    # Run Frank-Wolfe algorithms and alternating projectiopng, then record trajectory data
    println("\n\n\n ----------> Cyclic Block-coordinate vanilla CG")
    _, _, _, _, td_cyc_bc_cg = run_FW(config, FrankWolfe.CyclicUpdate(), prod_lmo)
    println("\n\n\n ----------> Cyclic Block-coordinate BPCG")
    _, _, _, _, td_cyc_bc_bpcg = run_FW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
    println("\n\n\n ----------> Full Block-coordinate BPCG")
    _, _, _, _, td_full_bc_cg = run_FW(config, FrankWolfe.FullUpdate(), FrankWolfe.BPCGStep(), prod_lmo)  
    println("\n\n\n ----------> Full BPCG")
    _, _, _, _, td_bpcg = run_FW(config, prod_lmo)    
    # println("\n\n\n ----------> AP")
    # _, _, _, _, td_ap = run_FW(config, prod_lmo, true)

    push_to_trajectories!(ni_flag, td_cyc_bc_cg, trajectories, primal)
    push_to_trajectories!(ni_flag, td_cyc_bc_bpcg, trajectories, primal)
    push_to_trajectories!(ni_flag, td_full_bc_cg, trajectories, primal)
    push_to_trajectories!(ni_flag, td_bpcg, trajectories, primal)
    # push_to_trajectories!(ni_flag, td_ap, trajectories, primal)

    # Save trajectories
    save_trajectories("examples/traj_$basename.jld2", trajectories_ni, trajectories_i)

    return trajectories
end

# Use parameters from YAML file
config = Config("examples/config.yml")
print_config(config)

# S = {x ∈ ℝⁿ₊, ∑x ≤ right_side}
lmo_unitsmplx = FrankWolfe.UnitSimplexOracle(5.0)
# S = {x ∈ ℝⁿ₊, ∑x = right_side}
lmo_probsmplx_1 = FrankWolfe.ProbabilitySimplexOracle(5.0)
lmo_probsmplx_2 = FrankWolfe.ProbabilitySimplexOracle(6.0)
bounds = generate_rand_float_vector(config)
# Polytope similar to ℓ-∞ ball with shifted bounds, centered in midpoint of [lowerbound, upperbound] for each dimension
lmo_infnormball_1 = FrankWolfe.ScaledBoundLInfNormBall(-ones(config.n), ones(config.n))
lmo_infnormball_2 = FrankWolfe.ScaledBoundLInfNormBall(5*ones(config.n), 7*ones(config.n))
lmo_infnormball_3 = FrankWolfe.ScaledBoundLInfNormBall(8*ones(config.n), 9*ones(config.n))
lmo_infnormball_4 = FrankWolfe.ScaledBoundLInfNormBall(4*ones(config.n), 9*ones(config.n))
# Polytope similar to ℓ₁ ball with shifted bounds, centered in midpoint of [lowerbound, upperbound] for each dimension
lmo_onenormball_1 = FrankWolfe.ScaledBoundL1NormBall(-ones(config.n), ones(config.n))
lmo_onenormball_2 = FrankWolfe.ScaledBoundL1NormBall(5*ones(config.n), 7*ones(config.n))
lmo_onenormball_3 = FrankWolfe.ScaledBoundL1NormBall(6*ones(config.n), 7*ones(config.n))
lmo_onenormball_4 = FrankWolfe.ScaledBoundL1NormBall(7*ones(config.n), 8*ones(config.n))

# Setup Linear Minimization Oracles for the polytopes
lmo_list_01 = [lmo_probsmplx_1, lmo_unitsmplx]     # intersecting
lmo_list_02 = [lmo_probsmplx_2, lmo_unitsmplx]     # non intersecting

lmo_list_03 = [lmo_unitsmplx, lmo_onenormball_1]   # intersecting
lmo_list_04 = [lmo_unitsmplx, lmo_onenormball_2]   # intersecting
lmo_list_05 = [lmo_unitsmplx, lmo_onenormball_3]   # non intersecting

lmo_list_06 = [lmo_unitsmplx, lmo_infnormball_1]   # intersecting
lmo_list_07 = [lmo_unitsmplx, lmo_infnormball_2]   # intersecting
lmo_list_08 = [lmo_unitsmplx, lmo_infnormball_4]   # intersecting
lmo_list_09 = [lmo_unitsmplx, lmo_infnormball_3]   # non intersecting

lmo_list_10 = [lmo_infnormball_2, lmo_onenormball_4]   # intersecting
lmo_list_11 = [lmo_infnormball_2, lmo_onenormball_3]   # intersecting
lmo_list_12 = [lmo_infnormball_1, lmo_onenormball_2]   # non intersecting

lmo_list = lmo_list_12
name, ni_flag = lmo_name(lmo_list)

# execute main
trajectories = main(config, lmo_list, ni_flag, "$(name)_k$(config.k)_n$(config.n)")

# Labels for the plots
labels = ["C-BC-FW", "C-BC-BPCG", "F-BC-BPCG", "F-BPCG"]#, "AP"]
# Plot trajectories
plot_trajectories(trajectories, labels, yscalelog=false, xscalelog=true, filename="examples/plot_$(name)_k$(config.k)_n$(config.n).png")