# `fw_algorithms.jl`
# Run Cyclic Block-Coordinate vanilla FW over product LMO
function run_FW(order, prod_lmo, x0, config)
    
    f = x -> objective(config, x)
    grad! = (storage, x) -> gradient!(config, storage, x)

    trajectories = []
    # x, v, primal_gap, dual_gap
    _, _, _, _, trajectory_data = FrankWolfe.block_coordinate_frank_wolfe(
        f,
        grad!,
        prod_lmo,
        x0,
        update_order=order,
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=config.max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
# Run Block-Coordinate BPCG with specific update order (full, cyclic, etc.) over product LMO
function run_FW(order, update_step, prod_lmo, x0, config)
    
    f = x -> objective(config, x)
    grad! = (storage, x) -> gradient!(config, storage, x)
    
    trajectories = []
    # x, v, primal_gap, dual_gap
    _, _, _, _, trajectory_data = FrankWolfe.block_coordinate_frank_wolfe(
        f,
        grad!,
        prod_lmo,
        x0,
        update_order=order,
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=config.max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        update_step=update_step,
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
# Run BPCG over full product LMO
function run_FW(prod_lmo, x0, config)
    
    f = x -> objective(config, x)
    grad! = (storage, x) -> gradient!(config, storage, x)
    
    trajectories = []
    # x, v, primal_gap, dual_gap
    _, _, _, _, _, trajectory_data = FrankWolfe.blended_pairwise_conditional_gradient(
        f,
        grad!,
        prod_lmo,
        x0,
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=config.max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
# Run Alternating Projections over product LMO
function run_FW(prod_lmo, x0, config, ap_flag)
    if ap_flag
        trajectories = []
        # x, v, dual_gap, infeasible
        _, _, _, _, trajectory_data = FrankWolfe.alternating_projections(
            prod_lmo,
            x0,
            epsilon=config.target_tolerance,
            max_iteration=config.max_iterations,
            memory_mode=FrankWolfe.InplaceEmphasis(),
            verbose=true,
            trajectory=true,
            print_iter=config.max_iterations / 10
        );
        push!(trajectories, trajectory_data);
    
        return trajectories
    end
end