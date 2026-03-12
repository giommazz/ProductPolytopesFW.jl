# `product_algorithms.jl`
using FrankWolfe
using LinearAlgebra

"""
    AwayStep <: FrankWolfe.UpdateStep

Update step for block-coordinate away-step Frank-Wolfe.
Stores lazy-step settings and an optional active set used by update rule.
"""
mutable struct AwayStep <: FrankWolfe.UpdateStep
    lazy::Bool
    active_set::Union{FrankWolfe.AbstractActiveSet, Nothing}
    lazy_tolerance::Float64
    # used to decide if moving towards FW vertex or away vertex in the lazy case
    phi::Float64
end

"""
    Base.copy(obj::AwayStep)

Copy an `AwayStep` instance, cloning the active set when present.
"""
function Base.copy(obj::AwayStep)
    if obj.active_set === nothing
        return AwayStep(obj.lazy, nothing, obj.lazy_tolerance, obj.phi)
    else
        return AwayStep(obj.lazy, copy(obj.active_set), obj.lazy_tolerance, obj.phi)
    end
end

# Constructors
AwayStep() = AwayStep(false, nothing, 2.0, Inf)
AwayStep(lazy::Bool) = AwayStep(lazy, nothing, 2.0, Inf)
const AWAYSTEP_RENORM_INTERVAL = 1000


"""
    FrankWolfe.update_block_iterate(...)

Block update rule for [`AwayStep`](@ref) in the block-coordinate setting.
Implements both lazy and non-lazy away/FW decisions and updates active set.
"""
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
        # Every `AWAYSTEP_RENORM_INTERVAL` iterations: remove tiny-weight atoms, then renormalize.
        renorm = mod(t, AWAYSTEP_RENORM_INTERVAL) == 0
        if away_step_taken # active set update when away step
            # `renorm` always true when `away_step_take` == true
            FrankWolfe.active_set_update!(s.active_set, -gamma, vertex_taken, true, index) # `vertex_taken` is `a`
        else # active set update when FW step
            FrankWolfe.active_set_update!(s.active_set, gamma, vertex_taken, renorm, index) # `vertex_taken` is `v`
        end
    end
    
    if mod(t, AWAYSTEP_RENORM_INTERVAL) == 0
        FrankWolfe.active_set_renormalize!(s.active_set)
    end

    x = FrankWolfe.muladd_memory_mode(memory_mode, x, gamma, d)
    return (dual_gap, vertex_taken, d, gamma, step_type)

end
 




"""
    run_BlockCoordinateFW(config, order, update_step, prod_lmo)

Run block-coordinate Frank-Wolfe on a product LMO with user-specified block
update order and update step (for example FW, BPCG, or AwayStep).
"""
function run_BlockCoordinateFW(
    config::Config,
    order::FrankWolfe.BlockCoordinateUpdateOrder, # FrankWolfe.CyclicUpdate(), FrankWolfe.FullUpdate()
    update_step::FrankWolfe.UpdateStep, # FrankWolfe.BPCGStep(), FrankWolfe.FrankWolfeStep(), AwayStep()
    prod_lmo::FrankWolfe.ProductLMO
    )

    # L-smoothness constant (used only when the chosen stepsize strategy is Shortstep).
    L = 1
    # Use the actual number of blocks in `prod_lmo`; some auxiliary calls may use fewer blocks than `config.k`.
    x0 = find_starting_point(config, prod_lmo)
    n_blocks = length(prod_lmo.lmos) 

    line_search = get_stepsize_strategy(config.stepsize_strategy, L)
    
    # Create one update-step instance per block.
    if update_step isa FrankWolfe.UpdateStep
        update_step = Tuple(copy(update_step) for _ in 1:n_blocks)
    end

    for (i, s) in enumerate(update_step)
        # `AwayStep` needs a block-local active set, `FrankWolfeStep` does not.
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

"""
    run_FullFW(config, FW_algorithm, prod_lmo)

Run a full-update Frank-Wolfe-style algorithm over the full product LMO.
The concrete algorithm is provided by `FW_algorithm`.
"""
function run_FullFW(
    config::Config,
    FW_algorithm::Function,
    prod_lmo::FrankWolfe.ProductLMO
    )

    # L-smoothness constant used by `get_stepsize_strategy` only when `stepsize_strategy == 1` (Shortstep).
    # set to 1 because our problem's L is 1
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

"""
    run_FullAFW(config, prod_lmo; memory_mode=FrankWolfe.InplaceEmphasis())

Run full away-step Frank-Wolfe (AFW) over the full product LMO.
Keeps a local trajectory callback to exclude post-processing states
(`ST_LAST` / `ST_POSTPROCESS`) from plotting logs.
"""
function run_FullAFW(
    config::Config,
    prod_lmo::FrankWolfe.ProductLMO;
    memory_mode::FrankWolfe.MemoryEmphasis=FrankWolfe.InplaceEmphasis(),
)
    # L-smoothness constant used by `get_stepsize_strategy` only when `stepsize_strategy == 1` (Shortstep).
    L = 1
    x0 = find_starting_point(config, prod_lmo)
    traj_data = Any[]
    internal_cb = let traj_data=traj_data
        function (state, active_set)
            if state.step_type != FrankWolfe.ST_LAST && state.step_type != FrankWolfe.ST_POSTPROCESS
                push!(
                    traj_data,
                    (state.t, state.primal, state.primal - state.dual_gap, state.dual_gap, state.time),
                )
            end
            return true
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
        trajectory=false, # we maintain `traj_data` ourselves via `internal_cb`
        callback=internal_cb,
    )

    if !isempty(traj_data)
        last_row = traj_data[end]
        return res.x, res.v, last_row[2], last_row[4], traj_data
    end
    return res.x, res.v, res.primal, res.dual_gap, traj_data
end

"""
    run_AlternatingProjections(config, prod_lmo, ap_flag)

Run alternating projections over a product LMO when `ap_flag` is `true`.
"""
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
