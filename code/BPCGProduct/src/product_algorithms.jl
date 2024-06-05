# `product_algorithms.jl`
# Run Cyclic Block-Coordinate vanilla FW over product LMO
function run_FW(config::Config, order::FrankWolfe.BlockCoordinateUpdateOrder, prod_lmo::FrankWolfe.ProductLMO)
    
    x0 = find_starting_point(config, prod_lmo)

    x, v, primal, fw_gap, trajectory_data = FrankWolfe.block_coordinate_frank_wolfe(
        objective,
        gradient!,
        prod_lmo,
        x0,
        update_order=order,
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=config.max_print_iterations,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );

    return x, v, primal, fw_gap, trajectory_data
end
# (Multiple dispatch) Run Block-Coordinate BPCG with specific update order (full, cyclic, etc.) over product LMO
function run_FW(config::Config, order::FrankWolfe.BlockCoordinateUpdateOrder, update_step::FrankWolfe.UpdateStep, prod_lmo::FrankWolfe.ProductLMO)
    
    x0 = find_starting_point(config, prod_lmo)

    x, v, primal, fw_gap, trajectory_data = FrankWolfe.block_coordinate_frank_wolfe(
        objective,
        gradient!,
        prod_lmo,
        x0,
        update_order=order,
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=config.max_print_iterations,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        update_step=update_step,
        verbose=true,
        trajectory=true,
    );  

    return x, v, primal, fw_gap, trajectory_data
end
# (Multiple dispatch) Run BPCG over full product LMO
function run_FW(config::Config, prod_lmo::FrankWolfe.ProductLMO)
    
    x0 = find_starting_point(config, prod_lmo)

    x, v, primal, fw_gap, trajectory_data = FrankWolfe.blended_pairwise_conditional_gradient(
        objective,
        gradient!,
        prod_lmo,
        x0,
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=config.max_print_iterations,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );

    return x, v, primal, fw_gap, trajectory_data
end
# (Multiple dispatch) Run Alternating Projections over product LMO
function run_FW(config::Config, prod_lmo::FrankWolfe.ProductLMO, ap_flag::Bool)
    
    if ap_flag
        
        # starting point is computed on only one set (e.g. the first LMO in the product)
        x0 = FrankWolfe.compute_extreme_point(prod_lmo.lmos[1], zeros(config.n))

        x, v, fw_gap, infeasible, trajectory_data = FrankWolfe.alternating_projections(
            prod_lmo,
            x0,
            epsilon=config.target_tolerance,
            max_iteration=config.max_iterations,
            #memory_mode=FrankWolfe.InplaceEmphasis(),
            verbose=true,
            trajectory=true,
            print_iter=config.max_print_iterations
        );
    
        return x, v, fw_gap, infeasible, trajectory_data
    end
end

# Push trajectory data into appropriate vector (intersecting or non-intersecting) based on the flag `ni_flag` 
function push_to_trajectories!(ni_flag::Bool, trajectory_data_curr::Vector{Any}, trajectories_ni::Vector{Any}, trajectories_i::Vector{Any}, primal::Float64)
    # `primal` ≠ 0.0: the polytopes don't intersect
    if ni_flag
        # Replace "Primal" with "Primal Gap" in the FW log, i.e., replace f(x) with f(x) - `primal` 
        trajectory_data_pg = compute_primal_gap(trajectory_data_curr, primal)
        push!(trajectories_ni, trajectory_data_pg)
    # `primal` == 0.0: the polytopes do intersect
    else    
        push!(trajectories_i, trajectory_data_curr)
    end
end
# (Multiple dispatch)
function push_to_trajectories!(ni_flag::Bool, trajectory_data_curr::Vector{Any}, trajectories::Vector{Any}, primal::Float64)
    # `primal` ≠ 0.0: the polytopes don't intersect
    if ni_flag
        # Replace "Primal" with "Primal Gap" in the FW log, i.e., replace f(x) with f(x) - `primal` 
        trajectory_data_pg = compute_primal_gap(trajectory_data_curr, primal)
        push!(trajectories, trajectory_data_pg)
    # `primal` == 0.0: the polytopes do intersect
    else    
        push!(trajectories, trajectory_data_curr)
    end
end

# Save trajectory data to given .jld2 file
function save_trajectories(filename::String, trajectories::Vector{Any})
    
    save(filename, Dict("trajectories" => trajectories))
    println("Saving data to $filename")
end
# (Multiple dispatch) handle several trajectory data elements
function save_trajectories(filename::String, trajectories...)
    
    dict = Dict{String, Vector{Any}}()
    for td in trajectories
        dict[@var_name(td)] = td
    end
    save(filename, dict)
    println("Saving data to $filename with multiple trajectory data")
end

# Load data from .jld2 file
function load_trajectories(filename::String)
    
    # Load the file
    f = load(filename) 
    # Extract the keys and values from the loaded dictionary
    trajectory_keys = keys(f)
    trajectory_values = [f[k] for k in trajectory_keys]
    
    println("Loaded data from $filename with keys: $keys")
    
    return trajectory_values
end