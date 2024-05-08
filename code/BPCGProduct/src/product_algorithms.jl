# `fw_algorithms.jl`
# Run BC (ALM setting) with vanilla FW
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
# Run BC (ALM setting) with BPCG
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
