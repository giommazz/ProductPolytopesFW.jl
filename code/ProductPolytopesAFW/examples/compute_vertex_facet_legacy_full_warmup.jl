# `compute_vertex_facet_legacy_full_warmup.jl`
# Run either within the Julia REPL as include("examples/compute_vertex_facet_legacy_full_warmup.jl")
# or from linux terminal with:
#   julia --project=. examples/compute_vertex_facet_legacy_full_warmup.jl > test.log 2>&1
#
# This script mirrors the structure of `compute_intersection_custom_full_warmup.jl` but builds
# instances directly from legacy LMOs to reproduce a vertex-facet intersection geometry.
using ProductPolytopesAFW
using FrankWolfe
using Plots

# ---------------------------------------------------------------------------------
# INSTANCE GENERATION (LEGACY LMOS)
# Build two 2-polytope instances in R^n:
#   - non-intersecting: box vs shifted diamond
#   - intersecting: box facet touched by one diamond vertex
#
# Geometry:
#   P1 = [0,1]^n (BoxLMO)
#   P2 = DiamondLMO with midpoint m and one vertex p on facet x1=1, others outside with x1>1
# This guarantees P1 ∩ P2 = {p}, with p in relint({x1=1} ∩ P1).
function build_vertex_facet_lmos(config::Config; alpha::Float64=0.2, beta::Float64=0.45, separation::Float64=0.5)
    if config.k != 2
        error("This script is defined for k=2. Found k=$(config.k).")
    end
    if config.n < 2
        error("This script requires n >= 2 to build a non-degenerate vertex-facet geometry.")
    end
    if !(0.0 < beta < 0.5)
        error("`beta` must satisfy 0 < beta < 0.5 so the touching vertex lies in the relative interior of the facet.")
    end
    if alpha <= 0.0
        error("`alpha` must be strictly positive.")
    end
    if separation <= 0.0
        error("`separation` must be strictly positive.")
    end

    # P1: unit box [0,1]^n
    lower_box = zeros(Float64, config.n)
    upper_box = ones(Float64, config.n)
    lmo_box = FrankWolfe.BoxLMO(lower_box, upper_box)

    # Intersecting P2 (diamond):
    # - midpoint m_j = 0.5 for j>=2 and m_1 = 1 + alpha
    # - vertex v1^- = (1, 0.5, ..., 0.5) touches relint of facet {x1=1}
    # - all other vertices have x1 > 1, hence lie outside the unit box
    lower_i = fill(0.5 - beta, config.n)
    upper_i = fill(0.5 + beta, config.n)
    lower_i[1] = 1.0
    upper_i[1] = 1.0 + 2.0 * alpha
    lmo_diamond_i = FrankWolfe.DiamondLMO(lower_i, upper_i)

    # Non-intersecting P2: same shape translated along e1
    lower_ni = copy(lower_i)
    upper_ni = copy(upper_i)
    lower_ni[1] += separation
    upper_ni[1] += separation
    lmo_diamond_ni = FrankWolfe.DiamondLMO(lower_ni, upper_ni)

    # Typed vectors to hit `create_product_lmo(::Vector{LinearMinimizationOracle})`
    lmos_ni = FrankWolfe.LinearMinimizationOracle[lmo_box, lmo_diamond_ni]
    lmos_i = FrankWolfe.LinearMinimizationOracle[lmo_box, lmo_diamond_i]

    # The touching point (expected unique intersection for the intersecting instance)
    p_touch = fill(0.5, config.n)
    p_touch[1] = 1.0

    return lmos_ni, lmos_i, p_touch
end

function compute_nonintersecting_opt(config::Config, lmos_ni::Vector{FrankWolfe.LinearMinimizationOracle})
    config_opt = modify_config(config, target_tolerance=config.target_tolerance_opt)
    _, _, primal_opt, _ = compute_distance(config_opt, lmos_ni)
    return primal_opt
end

# ---------------------------------------------------------------------------------
# MAIN FUNCTIONS
function repl_warmup(config::Config, lmo_pairs::Vector{Vector{FrankWolfe.LinearMinimizationOracle}})
    for (i, lmos) in enumerate(lmo_pairs)
        ni_flag = i == 1
        println()
        if ni_flag
            println("\nNon intersecting")
        else
            println("\nIntersecting")
        end
        prod_lmo = create_product_lmo(lmos)
        _, _, _, _, _ = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        _, _, _, _, _ = run_FullFW(config, FrankWolfe.frank_wolfe, prod_lmo)
        _, _, _, _, _ = run_FullFW(config, FrankWolfe.away_frank_wolfe, prod_lmo)
    end
end

function main(config::Config, lmo_pairs::Vector{Vector{FrankWolfe.LinearMinimizationOracle}}, opt_ni::Float64, labels::Vector{String})
    trajectories_ni, trajectories_i = [], []

    for (i, lmos) in enumerate(lmo_pairs)
        ni_flag = i == 1
        println()
        if ni_flag
            println("\nNon intersecting")
        else
            println("\nIntersecting")
        end

        prod_lmo = create_product_lmo(lmos)

        println("\n\n\n ----------> Cyclic Block-coordinate vanilla FW")
        _, _, _, _, td_cyc_bc_fw = run_BlockCoordinateFW(config, FrankWolfe.CyclicUpdate(), FrankWolfe.FrankWolfeStep(), prod_lmo)
        push_to_trajectories!(ni_flag, td_cyc_bc_fw, trajectories_ni, trajectories_i, opt_ni)

        println("\n\n\n ----------> Full FW")
        _, _, _, _, td_full_fw = run_FullFW(config, FrankWolfe.frank_wolfe, prod_lmo)
        push_to_trajectories!(ni_flag, td_full_fw, trajectories_ni, trajectories_i, opt_ni)

        println("\n\n\n ----------> Full AFW")
        _, _, _, _, td_full_afw = run_FullFW(config, FrankWolfe.away_frank_wolfe, prod_lmo)
        push_to_trajectories!(ni_flag, td_full_afw, trajectories_ni, trajectories_i, opt_ni)
    end

    return trajectories_ni, trajectories_i
end


t_start_tot = time()


# ---------------------------------------------------------------------------------
# YAML PARAMETERS
# ---------------------------------------------------------------------------------
println("Main configuration")
config = Config("examples/config.yml")
config = modify_config(config, k=2) # this script is defined for 2 blocks
print_config(config)

println("Warmup configuration")
warmup_n = max(2, min(config.n, 55))
config_warmup = modify_config(config, k=2, n=warmup_n)
print_config(config_warmup)


# ---------------------------------------------------------------------------------
# WARM-UP SCRIPT
# ---------------------------------------------------------------------------------
println()
println()
println("********************************************************")
println("WARMUP: Building legacy LMO instances")
println("********************************************************")
lmos_ni_warmup, lmos_i_warmup, p_touch_warmup = build_vertex_facet_lmos(config_warmup)
println("Warmup touching point p (first 5 coords): ", p_touch_warmup[1:min(5, length(p_touch_warmup))])

println("********************************************************")
println("WARMUP: Running FW on legacy instances")
println("********************************************************")
repl_warmup(config_warmup, [lmos_ni_warmup, lmos_i_warmup])


# ---------------------------------------------------------------------------------
# MAIN SCRIPT
# ---------------------------------------------------------------------------------
println()
println()
println()

results_dir = "examples/results_vertex_facet_legacy"
times_dir = ensure_dir(results_dir * "/times")
logs_dir = ensure_dir(results_dir * "/iter_logs")
plots_dir = ensure_dir(results_dir * "/plots")

labels = ["C-BC-FW", "F-FW", "F-AFW"]

println("********************************************************")
println("MAIN: Building legacy LMO instances and solving non-intersecting optimum")
println("********************************************************")
lmos_ni, lmos_i, p_touch = build_vertex_facet_lmos(config)
println("Touching point p (first 10 coords): ", p_touch[1:min(10, length(p_touch))])
opt_ni = compute_nonintersecting_opt(config, lmos_ni)
println("Non-intersecting optimal value (distance objective): ", opt_ni)

basename_run = "legacy_vertex_facet_" * generate_filename(config)

println("\n\n********************************************************")
println("MAIN: Running FW on legacy instances")
println("********************************************************")
t_start_main = time()
trajectories_ni, trajectories_i = main(config, [lmos_ni, lmos_i], opt_ni, labels)
t_end_main = time()
println("\n\n\t\tElapsed time main: ", t_end_main - t_start_main, " seconds\n\n")


# ---------------------------------------------------------------------------------
# PROCESS AND SAVE LOGS
# ---------------------------------------------------------------------------------
padded_trajectories_ni, min_length_ni, max_length_ni = pad_log_data(trajectories_ni)
padded_trajectories_i, min_length_i, max_length_i = pad_log_data(trajectories_i)

opt_ni = best_seen_solution(padded_trajectories_ni, opt_ni)

save_logdata_to_csv(padded_trajectories_ni, opt_ni, max_length_ni, labels, logs_dir, "ni_" * basename_run)
save_logdata_to_csv(padded_trajectories_i, max_length_i, labels, logs_dir, "i_" * basename_run)


# ---------------------------------------------------------------------------------
# CREATE AND SAVE PLOTS
# ---------------------------------------------------------------------------------
cutoff_trajectories_ni, _ = cutoff_log_shortest_time(padded_trajectories_ni)
cutoff_trajectories_i, _ = cutoff_log_shortest_time(padded_trajectories_i)

cutoff_trajectories_ni_pgap = [compute_primal_gap(t, opt_ni) for t in cutoff_trajectories_ni]

fig_ni = plot_time_only(config, cutoff_trajectories_ni_pgap, labels, yscalelog=true, xscalelog=true)
fig_i = plot_time_only(config, cutoff_trajectories_i, labels, yscalelog=true, xscalelog=true)

fig_ni_filename = plots_dir * "/plot_ni_$basename_run"
fig_i_filename = plots_dir * "/plot_i_$basename_run"

Plots.savefig(fig_ni, fig_ni_filename * ".pdf")
Plots.savefig(fig_i, fig_i_filename * ".pdf")


# ---------------------------------------------------------------------------------
# SAVE TIMES
# ---------------------------------------------------------------------------------
log_times(trajectories_i, labels, times_dir * "/times_i_" * basename_run)
log_times(trajectories_ni, labels, times_dir * "/times_ni_" * basename_run)


t_end_tot = time()
println("\n\n\t\tElapsed time tot: ", t_end_tot - t_start_tot, " seconds\n\n")