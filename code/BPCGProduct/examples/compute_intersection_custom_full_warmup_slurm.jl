# `compute_intersection_custom_instances.jl`
# Run either within the Julia REPL as include("/examples/compute_intersection_custom_full.jl")
# or from linux terminal with: `julia --project=. examples/compute_intersection_custom_full.jl > test.log 2>&1`
#
# This includes a little script that warms up the REPL: it is executed before running the main script, so that when the main script is run, 
#   the package `BPCGProduct` has already been precompiled and the plots obtained by `main` don't show initial arbitrary overhead
using BPCGProduct
using FrankWolfe
using Plots





# ---------------------------------------------------------------------------------
# MAIN FUNCTIONS
# Solve small instance (using `config_warmup`) to "warm-up" the REPL: this compiles `BPCGProduct`, so that no compilation needed upon running `main`
function repl_warmup(config::Config, vertices, shifted_vertices, primal, labels, basename)

    # Retrieve nonintersecting and intersecting LMOs from previously generated instances
    lmo_list = create_lmos(config, [vertices, shifted_vertices])

    # Will contain data about diafferent FW runs, for non-intersecting and intersecting polytopes
    trajectories_ni, trajectories_i = [], []

    for (i, lmos) in enumerate(lmo_list)
        # nonintersecting flag
        ni_flag = i == 1
        println()
        if ni_flag println("\nNon intersecting") else println("\nIntersecting") end
        prod_lmo = create_product_lmo(lmos)
        println("\n\n\n ----------> Full Block-coordinate Away FW (ours)")
        _, _, _, _, td_full_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), AwayStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_full_bc_afw, trajectories_ni, trajectories_i, primal)
    end

    return trajectories_ni, trajectories_i
end

function main(config::Config, vertices, shifted_vertices, primal, labels, basename)

    # Retrieve nonintersecting and intersecting LMOs from previously generated instances
    lmo_list = create_lmos(config, [vertices, shifted_vertices])

    # Will contain data about diafferent FW runs, for non-intersecting and intersecting polytopes
    trajectories_ni, trajectories_i = [], []

    for (i, lmos) in enumerate(lmo_list)
        # nonintersecting flag
        ni_flag = i == 1
        println()
        if ni_flag println("\nNon intersecting") else println("\nIntersecting") end

        prod_lmo = create_product_lmo(lmos)
        
        # Run Frank-Wolfe algorithms and alternating projections, then record trajectory data
        
        # ***************************************
        # Cyclic block-coordinate methods
        println("\n\n\n ----------> Cyclic Block-coordinate vanilla FW")
        _, _, _, _, td_cyc_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_cyc_bc_fw, trajectories_ni, trajectories_i, primal)
        
        # println("\n\n\n ----------> Cyclic Block-coordinate AFW")
        # _, _, _, _, td_cyc_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), AwayStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_cyc_bc_afw, trajectories_ni, trajectories_i, primal)

        # println("\n\n\n ----------> Stochastic Block-coordinate vanilla FW")
        # _, _, _, _, td_stoc_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.StochasticUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_stoc_bc_fw, trajectories_ni, trajectories_i, primal)
        
        # println("\n\n\n ----------> Stochastic Block-coordinate AFW")
        # _, _, _, _, td_stoc_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.StochasticUpdate(), AwayStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_stoc_bc_afw, trajectories_ni, trajectories_i, primal)
        
        # println("\n\n\n ----------> Cyclic Block-coordinate BPFW")
        # _, _, _, _, td_cyc_bc_bpfw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.BPCGStep(), prod_lmo)
        # push_to_trajectories!(ni_flag, td_cyc_bc_bpfw, trajectories_ni, trajectories_i, primal)

        # ***************************************
        # Full block-coordinate methods
        println("\n\n\n ----------> Full Block-coordinate vanilla FW")
        _, _, _, _, td_full_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_full_bc_fw, trajectories_ni, trajectories_i, primal)

        println("\n\n\n ----------> Full Block-coordinate Away FW (ours)")
        _, _, _, _, td_full_bc_afw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), AwayStep(), prod_lmo)  
        push_to_trajectories!(ni_flag, td_full_bc_afw, trajectories_ni, trajectories_i, primal)

        # println("\n\n\n ----------> Full Block-coordinate Blended Pairwise FW (ours)")
        # _, _, _, _, td_full_bc_bpfw = run_BlockCoordinateFW(config, FrankWolfe.FullUpdate(), FrankWolfe.BPCGStep(), prod_lmo)  
        # push_to_trajectories!(ni_flag, td_full_bc_bpfw, trajectories_ni, trajectories_i, primal)

        # ***************************************
        # Full methods
        println("\n\n\n ----------> Full FW")
        _, _, _, _, td_full_fw = run_FullFW(config, FrankWolfe.frank_wolfe, prod_lmo)    
        push_to_trajectories!(ni_flag, td_full_fw, trajectories_ni, trajectories_i, primal)

        println("\n\n\n ----------> Full AFW")
        _, _, _, _, td_full_afw = run_FullFW(config, FrankWolfe.away_frank_wolfe, prod_lmo)    
        push_to_trajectories!(ni_flag, td_full_afw, trajectories_ni, trajectories_i, primal)

        # println("\n\n\n ----------> Full BPFW")
        # _, _, _, _, td_full_bpfw = run_FullFW(config, FrankWolfe.blended_pairwise_conditional_gradient, prod_lmo)    
        # push_to_trajectories!(ni_flag, td_full_bpfw, trajectories_ni, trajectories_i, primal)
        
        # println("\n\n\n ----------> AP")
        # _, _, _, _, td_ap = run_AlternatingProjections(config, prod_lmo, true)    
        # # `FrankWolfe.alternating_projections` computes ||x-y|| rather than 1/2 ||x-y||
        # push_to_trajectories!(ni_flag, td_ap, trajectories_ni, trajectories_i, 2*primal)

        # Save trajectories
        # save_trajectories("examples/traj_$basename.jld2", trajectories_ni, trajectories_i)
    end

    return trajectories_ni, trajectories_i
end




# ---------------------------------------------------------------------------------
# YAML PARAMETERS
config = Config("examples/config.yml")
print_config(config)
config_warmup = modify_config(config, k=2, n=15)
print_config(config_warmup)
# ---------------------------------------------------------------------------------
# WARM-UP SCRIPT
# Generate instances 
println()
println()
println("********************************************************")
println("WARMUP: Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config_warmup)
# Optimal solution
primal = primal - 1     # Numerical reasons
basename = generate_filename(config_warmup)
# Labels for the plots
labels = ["F-BC-AFW"]
# execute main
println("********************************************************")
println("WARMUP: Running FW on the instances")
println("********************************************************")
_, _ = repl_warmup(config_warmup, vertices, shifted_vertices, primal, labels, basename)





# ---------------------------------------------------------------------------------
# MAIN SCRIPT
# ---------------------------------------------------------------------------------
println()
println()
println()
# ---------------------------------------------------------------------------------
# PARAMETERS
results_dir = "examples/results_linesearch_afw" # "results_shortstep", 
times_dir = results_dir*"/times"
logs_dir = results_dir*"/logs"
plots_dir = results_dir*"/plots"
# Labels for the plots
labels = ["C-BC-FW", "F-BC-FW", "F-BC-AFW", "F-FW", "F-AFW"] # ["C-BC-FW", "C-BC-AFW", "C-BC-BPFW", "F-BC-FW", "F-BC-AFW", "F-BC-BPFW", "F-FW", "F-AFW", "F-BPFW", "AP"]

# ---------------------------------------------------------------------------------
# GENERATE INSTANCES
println("********************************************************")
println("MAIN: Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config)
# Optimal solution
primal = primal - 1     # Numerical reasons
basename = generate_filename(config)

# ---------------------------------------------------------------------------------
# RUN MAIN
println("\n\n********************************************************")
println("MAIN: Running FW on the instances")
println("********************************************************")
trajectories_ni, trajectories_i = main(config, vertices, shifted_vertices, primal, labels, basename)


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
    trajectories::Vector{Any},
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
    # Custom colorblind palette (Okabe–Ito)
    okabe_ito = [
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
        color_choice = (length(trajectories) <= length(okabe_ito)) ?
                          okabe_ito[i] :
                          okabe_ito[mod1(i, length(okabe_ito))]
        
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

# Generate plots (do not pass `filename` argument, so .png is not automatically saved)
fig_ni = plot_trajectories(trajectories_ni, labels, yscalelog=true, xscalelog=true)
fig_i = plot_trajectories(trajectories_i, labels, yscalelog=true, xscalelog=true)
# fig_ni = plot_time_only(trajectories_ni, labels, yscalelog=true, xscalelog=true)
# fig_i  = plot_time_only(trajectories_i, labels, yscalelog=true, xscalelog=true)

# Decide filename
fig_ni_filename = plots_dir*"/plot_ni_$basename"
fig_i_filename = plots_dir*"/plot_i_$basename"
# Plot trajectories
# Plots.plot!(fig_ni, size=(1200, 800))  # Larger figure size
# Plots.plot!(fig_i, size=(1200, 800))  # Larger figure size
# Manually save only the PDF versions
Plots.savefig(fig_ni, fig_ni_filename*".pdf")  
Plots.savefig(fig_i, fig_i_filename*".pdf")

# ---------------------------------------------------------------------------------
# SAVE TIMES
# Save time data in `.csv` format
# log_data(trajectories_i, labels, times_dir*"/times_i_"*basename)
# log_data(trajectories_ni, labels, times_dir*"/times_ni_"*basename)


# println()
# for fw_variant_i in 1:length(trajectories_ni)
#     println("\t", trajectories_ni[fw_variant_i][1])
# end