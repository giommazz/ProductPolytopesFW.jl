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
# YAML PARAMETERS
config = Config("examples/config.yml")
print_config(config)
config_warmup = modify_config(config, k=2, n=15)
print_config(config_warmup)





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




println()
println()
# ---------------------------------------------------------------------------------
# WARM-UP SCRIPT
# Generate instances 
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




println()
println()
println()
# ---------------------------------------------------------------------------------
# MAIN SCRIPT
results_dir = "examples/results_linesearch_afw" # "results_shortstep", 
times_dir = results_dir*"/times"
logs_dir = results_dir*"/logs"
plots_dir = results_dir*"/plots"



# Generate instances 
println("********************************************************")
println("MAIN: Generating instances and solving them to optimum")
println("********************************************************")
vertices, shifted_vertices, primal, fw_gap = generate_polytopes(config)
# Optimal solution
primal = primal - 1     # Numerical reasons
basename = generate_filename(config)

# Labels for the plots
labels = ["C-BC-FW", "F-BC-FW", "F-BC-AFW", "F-FW", "F-AFW"] # ["C-BC-FW", "C-BC-AFW", "C-BC-BPFW", "F-BC-FW", "F-BC-AFW", "F-BC-BPFW", "F-FW", "F-AFW", "F-BPFW", "AP"]

# execute main
println("\n\n********************************************************")
println("MAIN: Running FW on the instances")
println("********************************************************")
trajectories_ni, trajectories_i = main(config, vertices, shifted_vertices, primal, labels, basename)


"""
function save_padded_logdata_to_csv(
    padded_trajectories::Vector{Vector{Tuple{T, T}}},
    max_length:Int64,
    labels::Vector{String},
    logs_dir::String,
    basename::String
    ) where T<:Number
    
    if length(labels) ≠ length(padded_trajectories)
        error("The number of labels ($(length(labels))) does not correspond to the number of FW algorithms run ($(length(padded_trajectories))).\nPlease fix this in your code.")
    end

    # transform labels by deleting dashes: `replace` takes a String and replaces the "-" symbol with ""
    labels_logs = [replace(l, "-" => "") for l in labels]

    # Store results in DataFrame object
    # todo: the dataframe will have one header with labels: iter, "labels_logs[fw_variant_i]*_pgap", "labels_logs[fw_variant_i]*_dual", "labels_logs[fw_variant_i]*_dgap", "labels_logs[fw_variant_i]*_time" for all of the FW variants, so for example if I ran 5 FW variants, there will be in total 1 + 4*5 labels
    times = DataFrames.DataFrame(iters=Int64[], ...)

    # each tuple/iteration records primal_gap (PG), dual (D), dual_gap (DG), time (T)
    for iter in 1:max_length
        row_curr_iter = ()
        for fw_variant_i in 1:length(padded_trajectories)
            todo: one row of the dataframe will contain, for each FW variant, 4 elements: 2-5 of the tuple trajectories_ni[fw_variant_i][iter]).
            so in general each dataframe row will have values corresponding to, for example: iter, cbcfw_pgap, cbcfw_dual, cbcfw_dgap, cbcfw_time, fbcfw_pgap, fbcfw_dual, ..., fafw_time
            fw_variant_subrow_curr_iter = ...
            row = push!(row_curr_iter, fw_variant_subrow_curr_iter)
        end
        # todo: push!(times, row_curr_iter)
    end
    # todo: CSV.write(logs_dir*"/"*basename*".csv", times)
end
"""

using DataFrames
using CSV


# Aggregate padded trajectories from multiple FW variants into a DataFrame and write it into a CSV file
function save_padded_logdata_to_csv(
    padded_trajectories::Vector{Vector{NTuple{5, T}}},
    max_length::Int64,
    labels::Vector{String},
    logs_dir::String,
    basename::String
    ) where T <: Number

    # Check that the number of labels matches the number of trajectories
    if length(labels) ≠ length(padded_trajectories)
        error("The number of labels ($(length(labels))) does not correspond to the number of FW algorithms run ($(length(padded_trajectories))).\nPlease fix this in your code.")
    end

    # Transform labels by deleting dashes.
    labels_logs = [replace(l, "-" => "") for l in labels]
    n_variants = length(padded_trajectories)

    # Prepare empty DataFrame. Labels: 'iter' + for each FW variant, 4 columns 'FWVar_pgap', 'FWVar_dual', 'FWVar_dgap', 'FWVar_time'
    df = DataFrame()
    
    # Create or initialize a column named `:iter` in the DataFrame
    # indexing with `!`, we access the column by reference ("inplace" update) → modify actual column rather than a copy
    # using `:` before `iter` defines it as a Symbol, ← column names in Julia DataFrames are typically Symbols
    # `Int64[]` builds the `:iter` column: an empty vector of type Int64
    df[!, :iter] = Int64[]
    
    # For every label in `labels_logs`, create 4 new columns
    for l in labels_logs
        # `Symbol(l * "_pgap")` concatenates `l` with "_pgap" and converts the result into a Symbol, which is then used as column name
        # `Vector{T}()` initializes an empty vector with element type T for that column
        df[!, Symbol(l * "_pgap")] = Vector{T}()
        df[!, Symbol(l * "_dual")] = Vector{T}()
        df[!, Symbol(l * "_dgap")] = Vector{T}()
        df[!, Symbol(l * "_time")] = Vector{T}()
    end

    # ------------------------------
    # Now we build the rows of the DataFrame.
    # Each row corresponds to one iteration (from 1 to max_length).
    # For each iteration, we extract a specific tuple from each FW variant.
    # Each tuple contains 5 elements, but we ignore the first element.
    # We take elements 2 to 5 (pgap, dual, dgap, time) and add them in order.
    # ------------------------------
    for iter in 1:max_length
        # Initialize a dictionary that will temporarily hold the data for the current row.
        # Keys will be Symbols corresponding to DataFrame column names, and values the data to put in that cell.
        row_dict = Dict{Symbol, Any}()
        row_dict[:iter] = iter
        # For each FW variant, extract its tuple for the current iteration index.
        for fw_variant_i in 1:n_variants
            # Extract the 5-tuple; expected format: (ignored, pgap, dual, dgap, time)
            t = padded_trajectories[fw_variant_i][iter]
            base = labels_logs[fw_variant_i]
            # Build the keys for each of the four associated columns.
            key_pgap = Symbol(base * "_pgap")
            key_dual = Symbol(base * "_dual")
            key_dgap = Symbol(base * "_dgap")
            key_time = Symbol(base * "_time")
            
            # Assign the corresponding tuple elements to the proper keys in the row dictionary.
            row_dict[key_pgap] = t[2]
            row_dict[key_dual] = t[3]
            row_dict[key_dgap] = t[4]
            row_dict[key_time] = t[5]
        end
        # ------------------------------
        # The following line pushes the current row (built as a NamedTuple) into the DataFrame.
        # Explanation:
        # - push!(df, ...): Adds a new row to the DataFrame.
        # - (; row_dict...): The semicolon here is used to create a NamedTuple from the dictionary.
        #   - The syntax (; row_dict...) uses the "splat" operator (three dots) to unpack the contents of row_dict as keyword arguments
        #     into the NamedTuple constructor. This is a concise way to convert a dictionary of key=>value pairs into a NamedTuple.
        #   The NamedTuple keys must match the DataFrame’s columns, which we ensured in our construction.
        # ------------------------------
        push!(df, (; row_dict...))
    end

    # Write the DataFrame to a CSV file.
    CSV.write(joinpath(logs_dir, basename * ".csv"), df)
end


# Save log data in `.csv` format
# pad data so that all FW runs have the same number of iterations/lines
padded_trajectories_ni, min_length_ni, max_length_ni = pad_log_data(trajectories_ni)
padded_trajectories_i, min_length_i, max_length_i = pad_log_data(trajectories_i)
readline()

# Save time data in `.csv` format
log_data(trajectories_i, labels, times_dir*"/times_i_"*basename)
log_data(trajectories_ni, labels, times_dir*"/times_ni_"*basename)

# Save plot trajectories in `.pdf` format
fig_ni_filename = plots_dir*"/plot_ni_$basename"
fig_i_filename = plots_dir*"/plot_i_$basename"
# Generate plots (do not pass `filename` argument, so .png is not automatically saved)
fig_ni = plot_trajectories(trajectories_ni, labels, yscalelog=true, xscalelog=true)
fig_i = plot_trajectories(trajectories_i, labels, yscalelog=true, xscalelog=true)
# Manually save only the PDF versions
#Plots.plot!(fig_ni, size=(1200, 800))  # Larger figure size
#Plots.plot!(fig_i, size=(1200, 800))  # Larger figure size
Plots.savefig(fig_ni, fig_ni_filename*".pdf")  
Plots.savefig(fig_i, fig_i_filename*".pdf")





println()
for fw_variant_i in 1:length(trajectories_ni)
    println("\t", trajectories_ni[fw_variant_i][1])
end
readline()