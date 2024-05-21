# `fw_algorithms.jl`
# Run Cyclic Block-Coordinate vanilla FW over product LMO
function run_FW(config::Config, order::FrankWolfe.BlockCoordinateUpdateOrder, prod_lmo::FrankWolfe.ProductLMO)
    
    f = x -> objective(config, x)
    grad! = (storage, x) -> gradient!(config, storage, x)
    x0 = find_starting_point(config, prod_lmo)
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
        print_iter=100, # config.max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
# (Multiple dispatch) Run Block-Coordinate BPCG with specific update order (full, cyclic, etc.) over product LMO
function run_FW(config::Config, order::FrankWolfe.BlockCoordinateUpdateOrder, update_step::FrankWolfe.UpdateStep, prod_lmo::FrankWolfe.ProductLMO)
    
    f = x -> objective(config, x)
    grad! = (storage, x) -> gradient!(config, storage, x)
    x0 = find_starting_point(config, prod_lmo)
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
        print_iter=100, # config.max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        update_step=update_step,
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
# (Multiple dispatch) Run BPCG over full product LMO
function run_FW(config::Config, prod_lmo::FrankWolfe.ProductLMO)
    
    f = x -> objective(config, x)
    grad! = (storage, x) -> gradient!(config, storage, x)
    x0 = find_starting_point(config, prod_lmo)
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
        print_iter=100, # config.max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
# TODO. DEBUG!!!
# (Multiple dispatch) Run Alternating Projections over product LMO
function run_FW(config::Config, prod_lmo::FrankWolfe.ProductLMO, ap_flag::Bool)
    if ap_flag
        trajectories = []
        x0 = find_starting_point(config, prod_lmo)

        # x, v, dual_gap, infeasible
        _, _, _, _, trajectory_data = FrankWolfe.alternating_projections(
            prod_lmo,
            x0,
            epsilon=config.target_tolerance,
            max_iteration=config.max_iterations,
            memory_mode=FrankWolfe.InplaceEmphasis(),
            verbose=true,
            trajectory=true,
            print_iter=100, # config.max_iterations / 10
        );
        push!(trajectories, trajectory_data);
    
        return trajectories
    end
end