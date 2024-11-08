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

    d = similar(x) # initialize direction vector

    # ********************
    # Step 0: compute away vertex, "local" FW vertex and related info (`lambda`, `loc`: active set weight and index)
    _, local_v, local_v_loc, _, a_lambda, a, a_loc, _, _ = FrankWolfe.active_set_argminmax(s.active_set, gradient) # away and "local" FW vertices
    grad_dot_x = fast_dot(gradient, x) # <∇f(xₜ), xₜ>
    grad_dot_a = fast_dot(gradient, a) # <∇f(xₜ), aₜ>
    away_gap = grad_dot_a - grad_dot_x # <∇f(xₜ), aₜ - xₜ>

    # Lazy version
    if s.lazy 

        # ********************
        # Step 1: compute step info (away vertex, local FW vertex, FW vertex)
        local_gap = grad_dot_x - grad_dot_local_v
        grad_dot_local_v = fast_dot(gradient, local_v) # <∇f(xₜ), local_vₜ>
        
        # lazy FW step
        if local_gap >= away_gap && local_gap >= s.phi / s.lazy_tolerance && local_gap >= epsilon    
            step_type = ST_LAZY
            gamma_max = one(a_lambda)
            d = FrankWolfe.muladd_memory_mode(memory_mode, d, x, local_v)
            vertex_taken = local_v
            away_step_taken = false
            fw_step_taken = true
            index = local_v_loc
        else
            # lazy away step
            if away_gap > local_gap && away_gap >= s.phi / s.lazy_tolerance
                step_type = ST_AWAY
                gamma_max = a_lambda / (1 - a_lambda)
                d = FrankWolfe.muladd_memory_mode(memory_mode, d, a, x)
                vertex_taken = a
                away_step_taken = true
                fw_step_taken = false
                index = a_loc
            # FW step (resort to calling the LMO)
            else
                
                # DEBUG: DA QUI...
                # optionally: try vertex storage
                if use_extra_vertex_storage
                    lazy_threshold = grad_dot_x - s.phi / s.lazy_tolerance
                    (found_better_vertex, new_forward_vertex) =
                        FrankWolfe.storage_find_argmin_vertex(extra_vertex_storage, gradient, lazy_threshold)
                    if found_better_vertex
                        @debug("Found acceptable lazy vertex in storage")
                        v = new_forward_vertex
                        step_type = ST_LAZYSTORAGE
                    else
                        v = FrankWolfe.compute_extreme_point(lmo, gradient)
                        step_type = ST_REGULAR
                    end
                else
                    v = FrankWolfe.compute_extreme_point(lmo, gradient)
                    step_type = ST_REGULAR
                end
                # DEBUG: ...A QUI

                # Real dual gap promises enough progress
                grad_dot_fw_vertex = fast_dot(v, gradient)
                dual_gap = grad_dot_x - grad_dot_fw_vertex
                if dual_gap >= s.phi / s.lazy_tolerance
                    gamma_max = one(a_lambda)
                    d = FrankWolfe.muladd_memory_mode(memory_mode, d, x, v)
                    fw_step_taken = true
                # Lower expectation for progress
                else
                    step_type = ST_DUALSTEP
                    s.phi = min(dual_gap, s.phi / 2.0)
                    gamma_max = zero(a_lambda)
                    fw_step_taken = false
                end
                vertex_taken = v
                away_step_taken = false
                index = -1
            end
        end

    # Non-lazy version
    else
        # ********************
        # Step 1: compute step info (away vertex, FW vertex)
        v = FrankWolfe.compute_extreme_point(lmo, gradient) # FW vertex
        dual_gap = grad_dot_x - fast_dot(gradient, v) # <∇f(xₜ), xₜ - vₜ>

        # FW step
        if dual_gap >= away_gap && dual_gap >= epsilon
            step_type = ST_REGULAR # memory-saving settings
            d = FrankWolfe.muladd_memory_mode(memory_mode, d, x, v)
            gamma_max = one(a_lambda)
            vertex_taken = v
            away_step_taken = false
            fw_step_taken = true
            index = -1 # will be pushed in active set
        
        # away step
        elseif away_gap >= epsilon
            step_type = ST_AWAY
            d =  FrankWolfe.muladd_memory_mode(memory_mode, d, a, x)
            gamma_max = a_lambda / (1 - a_lambda)
            vertex_taken = a
            away_step_taken = true
            fw_step_taken = false
            index = a_loc # atom weight will be updated in active set
        
        # no step, algorithm termination
        else
            step_type = ST_AWAY
            gamma_max = zero(a_lambda)
            vertex_taken = a
            away_step_taken = false
            fw_step_taken = false
            index = a_loc
        end
    end
    
    # ********************
    # Step 2: perform step and update active set
    gamma = 0.0
    # if a step is taken (i.e., algorithm will go on at least another iteration)
    if fw_step_taken || away_step_taken

        # perform stepsize strategy
        gamma = FrankWolfe.perform_line_search(
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
        gamma = min(gamma_max, gamma) # decide stepsize
        step_type = gamma ≈ gamma_max ? ST_DROP : step_type

        # DEBUG: DA QUI...
        # Update active set
        # every tot iterations: remove atoms w/small weight from active set, then renormalize weights to sum up to 1
        renorm = mod(t, s.renorm_interval) == 0
        if away_step_taken # active set update when away step
            # `renorm` always true when `away_step_take` == true
            # `add_dropped_vertices` decides if vertices dropped from active set are saved somewhere else
            FrankWolfe.active_set_update!(s.active_set, -gamma, vertex_taken, true, index) # `vertex_taken` is `a`
        else # active set update when FW step
            FrankWolfe.active_set_update!(s.active_set, gamma, vertex, renorm, index) # `vertex_taken` is `v`
        end
    end
    
    if mod(t, s.renorm_interval) == 0
        FrankWolfe.active_set_renormalize!(s.active_set)
        x = FrankWolfe.compute_active_set_iterate!(s.active_set)
    end

    x = FrankWolfe.muladd_memory_mode(memory_mode, x, gamma, d)
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