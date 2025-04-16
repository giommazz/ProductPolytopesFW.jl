# `plotting_utils.jl`

function compute_primal_gap(trajectories_curr::Vector{Any}, opt::Float64)
    trajectories_curr_pg = deepcopy(trajectories_curr)
    # Iterate over each tuple within the current block
    for i in 1:length(trajectories_curr_pg)
        # for FW algorithms:            iter, primal, dual,           dgap, time
        # for Alternating projections:  iter, infeas, partial infeas, dgap, time
        iter, primal, dual, dgap, time = trajectories_curr_pg[i]    
        # Compute primal gap
        pgap = primal - opt
        # Replace current tuple with a new one including the primal gap instead of the primal value
        trajectories_curr_pg[i] = (iter, pgap, dual, dgap, time)
    end
    return trajectories_curr_pg
end

"""
Given a vector of trajectories (each being a vector of 5‑tuples), extract time-based info and plot two subplots:
    - left: pgap (tuple field 2) over time (tuple field 5)
    - right: dgap (FW gap) (tuple field 4) over time

Arguments:
- `trajectories`: contains data to be plotted
- `labels`: as many as `length(trajectories)`
- `xscalelog`, `yscalelog`: axis scaling
- `primal_offset`: vertical offset for primal subplot (if needed).
- `offset`: starting index in each trajectory (default 2, as in the FrankWolfe.jl package)

Return: Plot object with two subplots
"""
function plot_time_only(
    trajectories::Vector{Vector{Any}},
    labels::Vector{String};
    xscalelog::Bool = true,   # Default to true for log scale on x-axis
    yscalelog::Bool = true,   # Default to true for log scale on y-axis
    legend_position = :topright,
    lstyle = fill(:solid, length(trajectories)),
    line_width = 1.3,
    primal_offset = 1e-8,
    offset = 2,
    size::Tuple{Int, Int} = (1200, 400)
)
    # Custom colorblind palette from https://www.color-hex.com/color-palette/49436
    colorblind_palette = [
        "#E69F00",  # Orange
        "#CC79A7",  # Pink
        "#0072B2",  # Blue
        "#009E73",  # Teal
        "#F0E442"   # Yellow
    ]

    # Create empty plot for the "primal gap" (pgap) over time.
    # We let Plots choose tick positions automatically (instead of hardcoding them)
    plt_primal = plot(
        xaxis = xscalelog ? :log10 : :identity,
        yaxis = yscalelog ? :log10 : :identity,
        xlabel = "Time (s)",
        ylabel = "Primal",
        legend = legend_position,
        lw = line_width,
        size = size,
        xguidefontsize = 12,
        yguidefontsize = 12,
    )

    # Create empty plot for the "dual gap" (FW gap) over time.
    plt_gap = plot(
        xaxis = xscalelog ? :log10 : :identity,
        yaxis = yscalelog ? :log10 : :identity,
        xlabel = "Time (s)",
        ylabel = "FW gap",
        legend = legend_position,
        lw = line_width,
        size = size
    )

    # Loop through each trajectory and add its data.
    # We check for sufficient length so that we don't pass an empty vector to the plotting routines.
    for (i, trajectory) in enumerate(trajectories)
        if length(trajectory) < offset
            continue  # Skip this trajectory if it does not contain enough points.
        end

        # Extract the data from the tuple assuming:
        # - Field 2 is the primal gap (pgap)
        # - Field 4 is the FW gap (dual gap)
        # - Field 5 is time.
        times = [trajectory[j][5] for j in offset:length(trajectory)]
        if isempty(times)
            continue  # Skip if no valid time points have been extracted.
        end
        primal_vals = [trajectory[j][2] + primal_offset for j in offset:length(trajectory)]
        gap_vals = [trajectory[j][4] for j in offset:length(trajectory)]
        
        # Choose a color (cycle through the palette if more series than colors are present)
        color_choice = (length(trajectories) <= length(colorblind_palette)) ?
                          colorblind_palette[i] :
                          colorblind_palette[mod1(i, length(colorblind_palette))]
        
        # Plot the current series on both subplots.
        plot!(plt_primal, times, primal_vals, label = labels[i],
              linestyle = lstyle[i], lw = line_width, color = color_choice)
        plot!(plt_gap, times, gap_vals, label = labels[i],
              linestyle = lstyle[i], lw = line_width, color = color_choice)
    end

    # Combine the two subplots side-by-side with a wide format
    final_plot = plot(plt_primal, plt_gap; layout = (1, 2), size = (1200, 400), bottom_margin = 8Plots.mm)
    return final_plot
end

# Input: each Vector in `trajectories` is a vector of 5-tuples: (iter, pgap, dual, dgap, time)
# Computed: `cutoff_time` = min of the final times among all runs
# The function returns new trajectories, truncated so that every run only contains entries where `time` ≤ `cutoff_time`
function cutoff_log_shortest_time(trajectories::Vector{Vector{Any}})
    
    # Decide `cutoff_time`: ∀ Vectors in `trajectories`, extract 5th element (time) of the last tuple, then compute min among all these times
    cutoff_time = minimum([last(traj)[5] for traj in trajectories])
    # Create `cutoff_trajectories`, truncated to earliest finish point: ∀ Vectors in `trajectories`, cut out tuples where `time` ≥ `cutoff_time`
    cutoff_trajectories = [
        filter(tuple -> tuple[5] ≤ cutoff_time, traj) for traj in trajectories
    ]
    
    return cutoff_trajectories, cutoff_time
end