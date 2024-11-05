# `product_algorithms.jl`

function compute_L(config::Config)
    return (config.k - 1) * sqrt(2 * config.k)
end

# Implementation Away-Step for block-coordinate setting
mutable struct AwayStep <: UpdateStep
    lazy::Bool
    active_set::Union{FrankWolfe.AbstractActiveSet, Nothing}
    lazy_tolerance::Float64
    # used to decide if moving towards FW vertex or away vertex in the lazy case
    phi::Float64
end

# Constructors
AwayStep() = AwayStep(false, nothing, 2.0)
AwayStep(lazy::Bool) = AwayStep(lazy, nothing, 2.0)

# Copy an AwayStep object with active set if it exists
function Base.copy(obj::AwayStep)
    if obj.active_set === nothing
        return AwayStep(nothing, obj.renorm_interval)
    else
        return AwayStep(copy(obj.active_set), obj.renorm_interval)
    end
end
# Update iterate function for AwayStep
function update_iterate(
    s::AwayStep,
    x,
    lmo,
    f,
    gradient,
    grad!,
    dual_gap,
    t,
    line_search,
    linesearch_workspace,
    memory_mode,
    epsilon,
    )

    d = similar(x) # initialize direction vector to zero

    # ********************
    # Step 1: active set FW vertex and away vertex
    # compute minimizer over active set, wrt gradient 
    # `a`, `a_lambda`, `a_loc`: away vertex, its active set weight and index
    _, _, _, _, a_lambda, a, a_loc, _, _ = FrankWolfe.active_set_argminmax(s.active_set, gradient)

    # <∇f(xₜ), aₜ>: scalar products between gradient and away vertex
    dot_away_vertex = fast_dot(gradient, a)
    # <∇f(xₜ), xₜ>: scalar products between gradient and current iterate
    grad_dot_x = fast_dot(gradient, x)
    # <∇f(xₜ), aₜ - xₜ>
    away_gap = dot_away_vertex - grad_dot_x

    # ********************
    # Step 2: FW vertex
    if !s.lazy
        v = FrankWolfe.compute_extreme_point(lmo, gradient) # compute FW vertex
        # <∇f(xₜ), xₜ - vₜ>
        dual_gap = grad_dot_x - fast_dot(gradient, v)
        s.phi = dual_gap
    end

    if dual_gap >= away_gap && dual_gap >= epsilon
        step_type = ST_REGULAR # memory-saving settings
        d = FrankWolfe.muladd_memory_mode(memory_mode, d, x, v)
        gamma_max = one(a_lambda)
        vertex_taken = v
        away_step_taken = false
        fw_step_taken = true
        index = -1
    
    elseif away_gap >= epsilon
        step_type = ST_AWAY
        d =  FrankWolfe.muladd_memory_mode(memory_mode, d, a, x)
        gamma_max = a_lambda / (1 - a_lambda)

        vertex_taken = a
        away_step_taken = true
        fw_step_taken = false
        index = a_loc
    else
        step_type = ST_AWAY
        gamma_max = zero(a_lambda)
        vertex_taken = a
        away_step_taken = false
        fw_step_taken = false
        index = a_loc
    end

    gamma = perform_line_search(
            line_search,
            t,
            f,
            grad!,
            gradient,
            x,
            d,
            gamma_max,
            linesearch_workspace,
            memory_mode,
        )
    gamma = min(gamma_max, gamma)
    
    step_type = gamma ≈ gamma_max ? ST_DROP : step_type

    # cleanup and renormalize every x iterations. Only for the FW steps
    renorm = mod(t, s.renorm_interval) == 0
    if away_step_taken
        FrankWolfe.active_set_update!(s.active_set, -gamma, vertex_taken, true, index, add_dropped_vertices=use_extra_vertex_storage, vertex_storage=extra_vertex_storage)
    else

        if add_dropped_vertices && gamma == gamma_max
            for vtx in s.active_set.atoms
                if vtx != v
                    push!(extra_vertex_storage, vtx)
                end
            end
        end
        active_set_update!(active_set, gamma, vertex, renorm, index)
    end

    x = muladd_memory_mode(memory_mode, x, gamma, d)
    return (dual_gap, vertex_taken, d, gamma, step_type)

end
 


# (Multiple dispatch) Run Block-Coordinate BPCG with specific update order (full, cyclic, etc.) over product LMO
function run_BlockCoordinateBlendedPairwiseFW(
    config::Config,
    order::FrankWolfe.BlockCoordinateUpdateOrder,
    update_step::FrankWolfe.UpdateStep, #FrankWolfe.BPCGStep(), FrankWolfe.FrankWolfeStep(), FrankWolfe.AwayStep()
    prod_lmo::FrankWolfe.ProductLMO
    )

    # L-smoothness constant
    L =  compute_L(config)

    x0 = find_starting_point(config, prod_lmo)

    x, v, primal, fw_gap, trajectory_data = FrankWolfe.block_coordinate_frank_wolfe(
        objective,
        gradient!,
        prod_lmo,
        x0,
        update_order=order,
        update_step=update_step,
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=get_stepsize_strategy(config.stepsize_strategy, L),
        print_iter=config.max_print_iterations,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );  

    return x, v, primal, fw_gap, trajectory_data
end
# (Multiple dispatch) Run BPCG over full product LMO
function run_FullBlendedPairwiseFW(
    config::Config,
    prod_lmo::FrankWolfe.ProductLMO
    )

    # L-smoothness constant
    L =  compute_L(config)

    x0 = find_starting_point(config, prod_lmo)

    x, v, primal, fw_gap, trajectory_data = FrankWolfe.blended_pairwise_conditional_gradient(
        objective,
        gradient!,
        prod_lmo,
        x0,
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=get_stepsize_strategy(config.stepsize_strategy, L),
        print_iter=config.max_print_iterations,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );

    return x, v, primal, fw_gap, trajectory_data
end
# (Multiple dispatch) Run Alternating Projections over product LMO
function run_AlternatingProjections(
    config::Config,
    prod_lmo::FrankWolfe.ProductLMO,
    ap_flag::Bool
    )
    
    if ap_flag
        
        # starting point is computed on only one set (e.g. the first LMO in the product)
        x0 = FrankWolfe.compute_extreme_point(prod_lmo.lmos[1], zeros(config.n))

        x, v, fw_gap, infeasible, trajectory_data = FrankWolfe.alternating_projections(
            prod_lmo,
            x0,
            epsilon=config.target_tolerance,
            max_iteration=config.max_iterations,
            memory_mode=FrankWolfe.InplaceEmphasis(),
            verbose=true,
            trajectory=true,
            print_iter=config.max_print_iterations
        );
    
        return x, v, fw_gap, infeasible, trajectory_data
    end
end

# Push trajectory data into appropriate vector (intersecting or non-intersecting) based on the flag `ni_flag` 
function push_to_trajectories!(
    ni_flag::Bool,
    trajectory_data_curr::Vector{Any},
    trajectories_ni::Vector{Any},
    trajectories_i::Vector{Any},
    primal::Float64
    )
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
function push_to_trajectories!(
    ni_flag::Bool,
    trajectory_data_curr::Vector{Any},
    trajectories::Vector{Any},
    primal::Float64
    )

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
function save_trajectories(
    filename::String,
    trajectories::Vector{Any}
    )
    
    save(filename, Dict("trajectories" => trajectories))
    println("Saving data to $filename")
end
# (Multiple dispatch) handle several trajectory data elements
function save_trajectories(
    filename::String,
    trajectories...
    )
    
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