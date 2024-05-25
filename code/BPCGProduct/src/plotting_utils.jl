# `plotting_utils.jl`

function compute_primal_gap(trajectories_curr::Vector{Any}, opt::Float64)
    trajectories_curr_pg = deepcopy(trajectories_curr)
    # Iterate over each tuple within the current block
    for i in 1:length(trajectories_curr_pg)
        iter, primal, dual, dgap, time = trajectories_curr_pg[i]    
        # Compute primal gap
        pgap = primal - opt
        # Replace current tuple with a new one including the primal gap instead of the primal value
        trajectories_curr_pg[i] = (iter, pgap, dual, dgap, time)
    end
    return trajectories_curr_pg
end
