#!/usr/bin/env julia
using CSV
using DataFrames
using Printf
using YAML

function usage()
    return "Usage: julia --project=. benchmarks_local_cvxh_lmo/summarize_results.jl [config.yml]"
end

function load_config(path::AbstractString)
    isfile(path) || error("Config file not found: $path")
    return YAML.load_file(path)
end

function results_dir(cfg_path::AbstractString, cfg)
    return abspath(joinpath(dirname(cfg_path), String(cfg["results_dir"])))
end

function parse_log_metrics(path::AbstractString)
    if !isfile(path)
        return (max_rss_kb=missing, exit_status=missing)
    end
    text = read(path, String)
    rss_match = match(r"Maximum resident set size \(kbytes\):\s+([0-9]+)", text)
    exit_match = match(r"Exit status:\s+([0-9]+)", text)
    rss = rss_match === nothing ? missing : parse(Int, rss_match.captures[1])
    exit_status = exit_match === nothing ? missing : parse(Int, exit_match.captures[1])
    return (max_rss_kb=rss, exit_status=exit_status)
end

function load_rows(csv_dir::AbstractString, raw_log_dir::AbstractString)
    rows = DataFrame()
    for csv_path in sort(filter(path -> endswith(path, ".csv"), readdir(csv_dir; join=true)))
        df = CSV.read(csv_path, DataFrame)
        nrow(df) == 1 || error("Expected one row in $(csv_path), found $(nrow(df)).")
        case_id = df.case_id[1]
        metrics = parse_log_metrics(joinpath(raw_log_dir, "$(case_id).log"))
        df.max_rss_kb = [metrics.max_rss_kb]
        df.exit_status = [metrics.exit_status]
        append!(rows, df; cols=:union)
    end
    return rows
end

function format_sci(value)
    if ismissing(value)
        return ""
    elseif value isa Real && isnan(value)
        return ""
    end
    return @sprintf("%.6e", Float64(value))
end

function compact_oracle_table(rows::DataFrame)
    selected = rows[rows.benchmark_kind .== "oracle", :]
    table = select(selected, :case_id, :n, :n_points, :iterations, :elapsed_s, :alloc_bytes, :alloc_count, :max_rss_kb)
    rename!(table, :iterations => :oracle_calls)
    table.elapsed_s = format_sci.(table.elapsed_s)
    table.alloc_bytes = format_sci.(table.alloc_bytes)
    table.alloc_count = format_sci.(table.alloc_count)
    table.max_rss_kb = format_sci.(table.max_rss_kb)
    sort!(table, [:n, :case_id])
    return table
end

function compact_algorithm_table(rows::DataFrame)
    selected = rows[in.(rows.benchmark_kind, Ref(["fw", "afw"])), :]
    table = select(selected, :case_id, :n, :n_points, :iterations, :elapsed_s, :alloc_bytes, :alloc_count, :max_rss_kb)
    table.elapsed_s = format_sci.(table.elapsed_s)
    table.alloc_bytes = format_sci.(table.alloc_bytes)
    table.alloc_count = format_sci.(table.alloc_count)
    table.max_rss_kb = format_sci.(table.max_rss_kb)
    sort!(table, [:n, :case_id])
    return table
end

function write_table_pair(base_path::AbstractString, table::DataFrame)
    csv_path = "$(base_path).csv"
    tsv_path = "$(base_path).tsv"
    CSV.write(csv_path, table)
    CSV.write(tsv_path, table; delim='\t')
    println("Wrote $(csv_path)")
    println("Wrote $(tsv_path)")
end

function main(args::Vector{String})
    cfg_path = if isempty(args)
        abspath(joinpath(@__DIR__, "config.yml"))
    elseif length(args) == 1
        abspath(args[1])
    else
        error(usage())
    end

    cfg = load_config(cfg_path)
    out_dir = results_dir(cfg_path, cfg)
    csv_dir = joinpath(out_dir, "csv")
    raw_log_dir = joinpath(out_dir, "raw_logs")
    isdir(csv_dir) || error("CSV directory not found: $csv_dir")

    rows = load_rows(csv_dir, raw_log_dir)
    sort!(rows, [:benchmark_kind, :n, :backend, :optimized, :cache])

    oracle_table = compact_oracle_table(rows)
    algorithm_table = compact_algorithm_table(rows)
    write_table_pair(joinpath(out_dir, "oracle_table"), oracle_table)
    write_table_pair(joinpath(out_dir, "algorithm_table"), algorithm_table)
end

try
    main(ARGS)
catch err
    showerror(stderr, err)
    println(stderr)
    println(stderr, usage())
    exit(1)
end
