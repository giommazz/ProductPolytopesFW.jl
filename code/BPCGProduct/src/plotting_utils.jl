# `plotting_utils.jl`

# Get, from each tuple: Iteration, Primal, Dual, Dual Gap, Time
function get_trajdata(iter::Tuple)
    return iter[1], iter[2], iter[3], iter[4], iter[5]
end

# Custom printing
function print_trajdata(trajdata::Vector{Any}, print_iter::Int64, opt::Float64)
    # `trajdata[1][1][i]` = (Iteration, Primal, Dual, Dual Gap, Time) → tuple w/data about i-th iteration of FW
    trajdata = trajdata[1][1]
    
    println("It      Primal              Primal Gap          Dual Gap            Time")
    println("-------------------------------------------------------------------------")
    
    # Print iterations
    for (idx, iter) in enumerate(trajdata)
        i, prim, _, dgap, time = get_trajdata(iter)
        pgap = max(0, prim - opt)  # Ensure pgap is not negative
        if idx == 1 || idx == length(trajdata) || idx % print_iter == 0
            @Printf.printf("%-8d %-18.6e %-18.6e %-18.6e %-18.6e\n", i, prim, pgap, dgap, time)
        end
    end
end

# return new trajectory data, replacing "Primal" f(x) with "Primal Gap" f(x)  -f(x*)
# function compute_primal_gap(trajectories_curr::Vector{Any}, opt::Float64)
#     trajectories_curr_pg = deepcopy(trajectories_curr)
#     # Iterate over each block of trajectory data
#     for i in eachindex(trajectories_curr_pg)
#         # Iterate over each tuple within the current block
#         for j in eachindex(trajectories_curr_pg[i])
#             iter, prim, dual, dgap, time = trajectories_curr_pg[i][j]
#             # Compute primal gap
#             pgap = prim - opt
#             # Replace current tuple with a new one including the primal gap instead of the primal value
#             trajectories_curr_pg[i][j] = (iter, pgap, dual, dgap, time)
#         end
#     end
#     return trajectories_curr_pg
# end


function compute_primal_gap(trajectories_curr::Vector{Any}, opt::Float64)
    trajectories_curr_pg = deepcopy(trajectories_curr)
    # Iterate over each block of trajectory data
    for i in eachindex(trajectories_curr_pg)
        # Iterate over each tuple within the current block
        for j in eachindex(trajectories_curr_pg[i])
            iter, prim, dual, dgap, time = trajectories_curr_pg[i][j]
            # Compute primal gap
            pgap = prim - opt
            # Replace current tuple with a new one including the primal gap instead of the primal value
            trajectories_curr_pg[i][j] = (iter, pgap, dual, dgap, time)
        end
    end
    return trajectories_curr_pg
end
