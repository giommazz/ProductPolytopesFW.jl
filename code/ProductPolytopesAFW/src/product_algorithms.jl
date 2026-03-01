# `product_algorithms.jl`
using FrankWolfe
using LinearAlgebra

# Implementation Away-Step for block-coordinate setting
mutable struct AwayStep <: FrankWolfe.UpdateStep
    lazy::Bool
    active_set::Union{FrankWolfe.AbstractActiveSet, Nothing}
    renorm_interval::Int
    lazy_tolerance::Float64
    # used to decide if moving towards FW vertex or away vertex in the lazy case
    phi::Float64
end

# Copy an AwayStep object with active set if it exists
function Base.copy(obj::AwayStep)
    if obj.active_set === nothing
        return AwayStep(obj.lazy, nothing, obj.renorm_interval, obj.lazy_tolerance, obj.phi)
    else
        return AwayStep(obj.lazy, copy(obj.active_set), obj.renorm_interval, obj.lazy_tolerance, obj.phi)
    end
end

# Constructors
AwayStep() = AwayStep(false, nothing, 1000, 2.0, Inf)
AwayStep(lazy::Bool) = AwayStep(lazy, nothing, 1000, 2.0, Inf)


# Update block iterate function for AwayStep
function FrankWolfe.update_block_iterate(
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
    d_container,
    )

    d = d_container # direction buffer

    # ********************
    # Step 0: compute away vertex, "local" FW vertex and related info (`lambda`, `loc`: active set weight and index)
    _, local_v, local_v_loc, _, a_lambda, a, a_loc, _, _ = FrankWolfe.active_set_argminmax(s.active_set, gradient) # away and "local" FW vertices
    grad_dot_x = dot(gradient, x) # <∇f(xₜ), xₜ>
    grad_dot_a = dot(gradient, a) # <∇f(xₜ), aₜ>
    away_gap = grad_dot_a - grad_dot_x # <∇f(xₜ), aₜ - xₜ>

    # Lazy version
    if s.lazy 

        # ********************
        # Step 1: compute step info (away vertex, local FW vertex, FW vertex)
        grad_dot_local_v = dot(gradient, local_v) # <∇f(xₜ), local_vₜ>
        local_gap = grad_dot_x - grad_dot_local_v
        
        # lazy FW step
        if local_gap >= away_gap && local_gap >= s.phi / s.lazy_tolerance && local_gap >= epsilon    
            step_type = FrankWolfe.ST_LAZY
            gamma_max = one(a_lambda)
            d = FrankWolfe.muladd_memory_mode(memory_mode, d, x, local_v)
            vertex_taken = local_v
            away_step_taken = false
            fw_step_taken = true
            index = local_v_loc
        else
            # lazy away step
            if away_gap > local_gap && away_gap >= s.phi / s.lazy_tolerance
                step_type = FrankWolfe.ST_AWAY
                gamma_max = a_lambda / (1 - a_lambda)
                d = FrankWolfe.muladd_memory_mode(memory_mode, d, a, x)
                vertex_taken = a
                away_step_taken = true
                fw_step_taken = false
                index = a_loc
            # FW step (resort to calling the LMO)
            else
                v = FrankWolfe.compute_extreme_point(lmo, gradient)
                step_type = FrankWolfe.ST_REGULAR

                # Real dual gap promises enough progress
                grad_dot_fw_vertex = dot(v, gradient)
                dual_gap = grad_dot_x - grad_dot_fw_vertex
                if dual_gap >= s.phi / s.lazy_tolerance
                    gamma_max = one(a_lambda)
                    d = FrankWolfe.muladd_memory_mode(memory_mode, d, x, v)
                    fw_step_taken = true
                # Lower expectation for progress
                else
                    step_type = FrankWolfe.ST_DUALSTEP
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
        dual_gap = grad_dot_x - dot(gradient, v) # <∇f(xₜ), xₜ - vₜ>

        # FW step
        if dual_gap >= away_gap && dual_gap >= epsilon
            step_type = FrankWolfe.ST_REGULAR # memory-saving settings
            d = FrankWolfe.muladd_memory_mode(memory_mode, d, x, v)
            gamma_max = one(a_lambda)
            vertex_taken = v
            away_step_taken = false
            fw_step_taken = true
            index = -1 # will be pushed in active set
        
        # away step
        elseif away_gap >= epsilon
            step_type = FrankWolfe.ST_AWAY
            d =  FrankWolfe.muladd_memory_mode(memory_mode, d, a, x)
            gamma_max = a_lambda / (1 - a_lambda)
            vertex_taken = a
            away_step_taken = true
            fw_step_taken = false
            index = a_loc # atom weight will be updated in active set
        
        # no step, algorithm termination
        else
            step_type = FrankWolfe.ST_AWAY
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
        step_type = gamma ≈ gamma_max ? FrankWolfe.ST_DROP : step_type

        # Update active set
        # Every `renorm_interval` iterations: remove atoms w/small weight from active set, then renormalize weights to sum up to 1
        renorm = mod(t, s.renorm_interval) == 0
        if away_step_taken # active set update when away step
            # `renorm` always true when `away_step_take` == true
            FrankWolfe.active_set_update!(s.active_set, -gamma, vertex_taken, true, index) # `vertex_taken` is `a`
        else # active set update when FW step
            FrankWolfe.active_set_update!(s.active_set, gamma, vertex_taken, renorm, index) # `vertex_taken` is `v`
        end
    end
    
    if mod(t, s.renorm_interval) == 0
        FrankWolfe.active_set_renormalize!(s.active_set)
    end

    x = FrankWolfe.muladd_memory_mode(memory_mode, x, gamma, d)
    return (dual_gap, vertex_taken, d, gamma, step_type)

end
 




# (Multiple dispatch) Run Block-Coordinate BPCG with specific update order (full, cyclic, etc.) over product LMO
function run_BlockCoordinateFW(
    config::Config,
    order::FrankWolfe.BlockCoordinateUpdateOrder, # FrankWolfe.CyclicUpdate(), FrankWolfe.FullUpdate()
    update_step::FrankWolfe.UpdateStep, # FrankWolfe.BPCGStep(), FrankWolfe.FrankWolfeStep(), AwayStep()
    prod_lmo::FrankWolfe.ProductLMO
    )

    # L-smoothness constant
    L = 1
    # DEBUG: notice that I couldn't use config.k because I sometimes call the function on two sets only
    x0 = find_starting_point(config, prod_lmo)
    n_blocks = length(prod_lmo.lmos) 

    line_search = get_stepsize_strategy(config.stepsize_strategy, L)
    
    # DEBUG: for some reason I had to actually convert the stuff below into a Tuple...which doesn't seem to be necessary in the block_coordinate_algorithms.jl in the package...
    if update_step isa FrankWolfe.UpdateStep
        update_step = Tuple(copy(update_step) for _ in 1:n_blocks)
    end
    if line_search isa FrankWolfe.LineSearchMethod
        line_search = Tuple(line_search for _ in 1:n_blocks)
    end

    for (i, s) in enumerate(update_step)
        # DEBUG: here I used `!(s isa FrankWolfeStep)` bc I am sure that vanillaFW doesn't use an active set. But I don't know which other updatesteps DON'T use an active set..
        if !(s isa FrankWolfe.FrankWolfeStep) && s.active_set === nothing
            s.active_set = FrankWolfe.ActiveSet([(1.0, copy(x0.blocks[i]))])
        end
    end
    
    res = FrankWolfe.block_coordinate_frank_wolfe(
        convex_feasibility_objective_v2b,
        convex_feasibility_gradient_v2!,
        prod_lmo,
        x0,
        update_order=order,
        update_step=update_step,
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=line_search,
        print_iter=config.max_print_iterations,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );  

    return res.x, res.v, res.primal, res.dual_gap, res.traj_data
end

# Run BPCG over full product LMO
function run_FullFW(
    config::Config,
    FW_algorithm::Function,
    prod_lmo::FrankWolfe.ProductLMO
    )

    # L-smoothness constant
    L = 1

    x0 = find_starting_point(config, prod_lmo)

    res = FW_algorithm(
        convex_feasibility_objective_v2b,
        convex_feasibility_gradient_v2!,
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

    return res.x, res.v, res.primal, res.dual_gap, res.traj_data
end

# Run full away-step Frank-Wolfe (AFW) over the full product LMO, exposing AFW-specific controls.
#
# Motivation:
# - `run_FullFW(..., FrankWolfe.away_frank_wolfe, ...)` hides AFW keyword arguments such as
#   `renorm_interval` and `weight_purge_threshold`, helpful for active-set management.
# - FrankWolfe.jl's default trajectory callback for AFW may include post-processing states. For debugging
#   we maintain our own trajectory and (by default) skip `ST_LAST` / `ST_POSTPROCESS`.
function run_FullAFW(
    config::Config,
    prod_lmo::FrankWolfe.ProductLMO;
    renorm_interval::Int=config.max_iterations + 10,
    weight_purge_threshold::Union{Nothing,Float64}=1e-32,
    memory_mode::FrankWolfe.MemoryEmphasis=FrankWolfe.InplaceEmphasis(),
    callback::Union{Nothing,Function}=nothing,
    debug_callback::Bool=false,
    spike_factor::Float64=10.0,
    log_postprocess::Bool=false,
    recompute_metrics::Bool=false,
)
    renorm_interval >= 1 || error("`renorm_interval` must be ≥ 1, got $renorm_interval")
    spike_factor > 1 || error("`spike_factor` must be > 1, got $spike_factor")

    # L-smoothness constant used by `get_stepsize_strategy` when `stepsize_strategy == 1`.
    L = 1

    x0 = find_starting_point(config, prod_lmo)

    if weight_purge_threshold === nothing
        weight_purge_threshold = FrankWolfe.weight_purge_threshold_default(eltype(x0))
    end

    # Custom trajectory container: store (t, primal, dual, dual_gap, time) tuples.
    traj_data = Any[]
    grad_buf = similar(x0)
    prev_primal = Ref{Float64}(Inf)

    internal_cb = let traj_data=traj_data,
        grad_buf=grad_buf,
        prev_primal=prev_primal,
        debug_callback=debug_callback,
        spike_factor=spike_factor,
        log_postprocess=log_postprocess,
        recompute_metrics=recompute_metrics

        function (state, active_set)
            # Keep post-processing states out of the main trajectory by default, but allow printing them.
            if state.step_type == FrankWolfe.ST_LAST || state.step_type == FrankWolfe.ST_POSTPROCESS
                if log_postprocess
                    println(
                        "AFW postprocess: t=$(state.t) step=$(state.step_type) primal=$(state.primal) dgap=$(state.dual_gap) |AS|=$(length(active_set))",
                    )
                end
                return true
            end

            # Recompute primal (always) at the current iterate `state.x`.
            x_curr = state.x
            primal_now = convex_feasibility_objective_v2b(x_curr)

            dual_gap_now = state.dual_gap
            if debug_callback
                w = active_set.weights
                wsum = sum(w)
                wmin = minimum(w)
                wmax = maximum(w)

                spike = isfinite(prev_primal[]) && primal_now > spike_factor * prev_primal[]
                invalid_as = !(wsum ≈ 1.0) || (wmin < -1e-12) || !isfinite(wsum) || !isfinite(wmin) || !isfinite(wmax)

                # On suspicious events (or if explicitly requested), recompute gradient and FW gap at `x_curr`
                # to keep logged primal/dual_gap consistent for diagnosis.
                if recompute_metrics || spike || invalid_as
                    convex_feasibility_gradient_v2!(grad_buf, x_curr)
                    v_now = FrankWolfe.compute_extreme_point(state.lmo, grad_buf)
                    dual_gap_now = dot(grad_buf, x_curr) - dot(grad_buf, v_now)
                end

                if spike || invalid_as
                    println(
                        "AFW DEBUG: t=$(state.t) step=$(state.step_type) γ=$(state.gamma) primal=$(primal_now) dgap=$(dual_gap_now) |AS|=$(length(active_set)) wsum=$(wsum) wmin=$(wmin) wmax=$(wmax)",
                    )
                end
                prev_primal[] = primal_now
            elseif recompute_metrics
                # If debugging is disabled but a full recompute was explicitly requested, do it.
                convex_feasibility_gradient_v2!(grad_buf, x_curr)
                v_now = FrankWolfe.compute_extreme_point(state.lmo, grad_buf)
                dual_gap_now = dot(grad_buf, x_curr) - dot(grad_buf, v_now)
            end

            push!(
                traj_data,
                (state.t, primal_now, primal_now - dual_gap_now, dual_gap_now, state.time),
            )

            return true
        end
    end

    # If the caller passed an extra callback (e.g., for debugging / early stopping),
    # run it after our internal trajectory logging callback.
    cb = if callback === nothing
        internal_cb
    else
        (state, active_set) -> begin
            cont = internal_cb(state, active_set)
            cont && callback(state, active_set)
        end
    end

    res = FrankWolfe.away_frank_wolfe(
        convex_feasibility_objective_v2b,
        convex_feasibility_gradient_v2!,
        prod_lmo,
        x0;
        epsilon=config.target_tolerance,
        max_iteration=config.max_iterations,
        line_search=get_stepsize_strategy(config.stepsize_strategy, L),
        print_iter=config.max_print_iterations,
        memory_mode=memory_mode,
        verbose=true,
        trajectory=false, # we maintain `traj_data` ourselves via `cb`
        callback=cb,
        renorm_interval=renorm_interval,
        weight_purge_threshold=weight_purge_threshold,
    )

    # Return last logged values as primal/gap if available (more consistent with the logged trajectory).
    if !isempty(traj_data)
        last_row = traj_data[end]
        return res.x, res.v, last_row[2], last_row[4], traj_data
    end
    return res.x, res.v, res.primal, res.dual_gap, traj_data
end

# Run Alternating Projections over product LMO
function run_AlternatingProjections(
    config::Config,
    prod_lmo::FrankWolfe.ProductLMO,
    ap_flag::Bool
    )
    
    if ap_flag
        
        # starting point is computed on only one set (e.g. the first LMO in the product)
        x0 = FrankWolfe.compute_extreme_point(prod_lmo.lmos[1], zeros(config.n))

        res = FrankWolfe.alternating_projections(
            prod_lmo,
            x0,
            epsilon=config.target_tolerance,
            max_iteration=config.max_iterations,
            memory_mode=FrankWolfe.InplaceEmphasis(),
            verbose=true,
            trajectory=true,
            print_iter=config.max_print_iterations
        );
    
        return res.x, res.v, res.dual_gap, res.status != FrankWolfe.STATUS_OPTIMAL, res.traj_data
    end
end
