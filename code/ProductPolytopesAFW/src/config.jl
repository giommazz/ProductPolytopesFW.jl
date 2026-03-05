# `config.jl`
struct Config
    
    # Number of polytopes
    k::Int
    # Polytope dimension
    n::Int
    # Either 0 or a list of `k` integers, each indicating the number of points to be used to generate intersecting polytopes
    n_points::Union{Int, Vector{Int}}
    # Upper bound multiplier used only when `n_points == 0` (auto-generate `n_points[i]` in {n+1, …, n*ub_n_points})
    ub_n_points::Int
    # Vertex sampling strategy used to generate polytope points/atoms (e.g., "box_uniform", "sphere", or "ellipsoid")
    vertex_sampling::String
    # Used only when `vertex_sampling == "sphere"` or `"ellipsoid"`: factor ρ ∈ (0, 1) controlling the inscribed sphere/ellipsoid size
    sphere_radius_factor::Float64
    # epsilon-optimality threshold
    target_tolerance::Float64
    # epsilon-optimality threshold for finding optimal solution
    target_tolerance_opt::Float64
    # Number of FW iterations
    max_iterations::Int
    # Number of FW iterations to compute optimal solution
    max_iterations_opt::Int
    # How often FW iteration log is printed to screen
    max_print_iterations::Int
    # Set the seed for reproducibility
    seed::Int
    # Toggle diagnostic prints (e.g., in `generate_polytopes`)
    verbose::Bool
    # Use convex-hull LMOs (true) or FrankWolfe.MathOptLMO (false)
    cvxhflag::Bool
    # Backend used when `cvxhflag == true`:
    # - "matrix": ProductPolytopesAFW.MatrixConvexHullLMO
    # - "vector": FrankWolfe.ConvexHullLMO (best effort: row views, no full copy)
    convex_hull_backend::String
    # MatrixConvexHullLMO cache cap:
    # -1 => auto cap, 0 => cache off, >0 => fixed cap
    matrix_lmo_cache_cap::Int
    # Use optimized extreme-point search (GEMV + argmin) inside MatrixConvexHullLMO
    matrix_lmo_use_optimized_search::Bool
    # Intersection anchor: point in Rⁿ where polytopes are aligned (see examples/config.yml for allowed values)
    intersection_anchor::String
    # Reference point used for each polytope when aligning to the anchor (e.g., "center" or "vertex")
    intersection_reference_point::String
    # Scalar parameter t used for certain anchor types (e.g., "p1_p2_segment"); non‑negative
    intersection_anchor_t::Float64
    # stepsize strategy:
    # - 0: SafeGoldenratio ("goldensearch", safer gamma reconstruction)
    # - 1: Shortstep with given L (specified in `product_algorithms.jl`)
    stepsize_strategy::Int 
end

# Constructor with default values
function Config()
    
    # Default intersection settings
    c = Config(
        2,              # k
        2,              # n
        [5, 5],         # n_points
        10,             # ub_n_points
        "box_uniform",  # vertex_sampling
        0.95,           # sphere_radius_factor
        1e-6,           # target_tolerance
        1e-08,          # target_tolerance_opt
        1000,           # max_iterations
        5000,           # max_iterations_opt
        100,            # max_print_iterations
        42,             # seed
        true,           # verbose
        true,           # cvxhflag
        "matrix",       # convex_hull_backend
        -1,             # matrix_lmo_cache_cap
        true,           # matrix_lmo_use_optimized_search
        "origin",       # intersection_anchor
        "center",       # intersection_reference_point
        0.0,            # intersection_anchor_t
        0,              # stepsize_strategy
    )
    
    return c 
end

# Constructor from YAML file and possibly set values
# Can be called, for eaxmple, as Config(<yaml_filename>) or Config(<yaml_filename>; n=n, k=k)
function Config(yaml_file::String; kwargs...)
    
    # Load and validate configuration from YAML file
    config = YAML.load(open(yaml_file))

    for (key, value) in kwargs
        config[string(key)] = value
    end

    # Validation checks (for YAML file dict)
    validate_config(config)

    # Handle the n_points field
    n_points = config["n_points"]
    ub_n_points = get(config, "ub_n_points", 10)
    if n_points == 0
        n_points = generate_n_points(config["n"], config["k"], config["seed"]; ub_n_points=ub_n_points)
    end
    vertex_sampling = get(config, "vertex_sampling", "box_uniform")
    sphere_radius_factor = get(config, "sphere_radius_factor", 0.95)
    verbose = get(config, "verbose", true)

    intersection_cfg = config["intersection"]
    anchor_t = get(intersection_cfg, "anchor_t", 0.0)

    c = Config(
        config["k"],
        config["n"],
        n_points,
        ub_n_points,
        vertex_sampling,
        sphere_radius_factor,
        config["target_tolerance"],
        config["target_tolerance_opt"],
        config["max_iterations"],
        config["max_iterations_opt"],
        config["max_print_iterations"],
        config["seed"],
        verbose,
        config["cvxhflag"],
        get(config, "convex_hull_backend", "matrix"),
        get(config, "matrix_lmo_cache_cap", -1),
        get(config, "matrix_lmo_use_optimized_search", true),
        intersection_cfg["anchor"],
        intersection_cfg["reference_point"],
        anchor_t,
        config["stepsize_strategy"],
    )
    return c 
end

# Function to update Config with any number of arguments
function modify_config(config::Config; kwargs...)
    
    # `get(kwargs, :key, default_value)` retrieves value for `:key` from `kwargs` if it exists, o/w returns default_value    
    k = get(kwargs, :k, config.k)
    n = get(kwargs, :n, config.n)
    n_points = get(kwargs, :n_points, config.n_points)
    ub_n_points = get(kwargs, :ub_n_points, config.ub_n_points)
    vertex_sampling = get(kwargs, :vertex_sampling, config.vertex_sampling)
    sphere_radius_factor = get(kwargs, :sphere_radius_factor, config.sphere_radius_factor)
    target_tolerance = get(kwargs, :target_tolerance, config.target_tolerance)
    target_tolerance_opt = get(kwargs, :target_tolerance_opt, config.target_tolerance_opt)
    max_iterations = get(kwargs, :max_iterations, config.max_iterations)
    max_iterations_opt = get(kwargs, :max_iterations_opt, config.max_iterations_opt)
    max_print_iterations = get(kwargs, :max_iterations, config.max_print_iterations)
    seed = get(kwargs, :seed, config.seed)
    verbose = get(kwargs, :verbose, config.verbose)
    cvxhflag = get(kwargs, :cvxhflag, config.cvxhflag)
    convex_hull_backend = get(kwargs, :convex_hull_backend, config.convex_hull_backend)
    matrix_lmo_cache_cap = get(kwargs, :matrix_lmo_cache_cap, config.matrix_lmo_cache_cap)
    matrix_lmo_use_optimized_search = get(kwargs, :matrix_lmo_use_optimized_search, config.matrix_lmo_use_optimized_search)
    intersection_anchor = get(kwargs, :intersection_anchor, config.intersection_anchor)
    intersection_reference_point = get(kwargs, :intersection_reference_point, config.intersection_reference_point)
    intersection_anchor_t = get(kwargs, :intersection_anchor_t, config.intersection_anchor_t)
    stepsize_strategy = get(kwargs, :stepsize_strategy, config.stepsize_strategy)

    # if `n_points` was supplied...
    if haskey(kwargs, :n_points)
        n_points = kwargs[:n_points] # what the caller passed
    else
        # ...else: fall-back, generate a fresh list that matches (k, n, seed)
        n_points = generate_n_points(n, k, seed; ub_n_points=ub_n_points)
    end

    # Return a new Config object with updated values
    c = Config(
        k,
        n,
        n_points,
        ub_n_points,
        vertex_sampling,
        sphere_radius_factor,
        target_tolerance,
        target_tolerance_opt,
        max_iterations,
        max_iterations_opt,
        max_print_iterations,
        seed,
        verbose,
        cvxhflag,
        convex_hull_backend,
        matrix_lmo_cache_cap,
        matrix_lmo_use_optimized_search,
        intersection_anchor,
        intersection_reference_point,
        intersection_anchor_t,
        stepsize_strategy,
    )
    
    return c
end

# Validation checks (for YAML file dict)
function validate_config(yaml_config::Dict{Any, Any})
    
    # Check for `k`
    if typeof(yaml_config["k"]) != Int || yaml_config["k"] < 1
        error("Invalid value for 'k': must be a positive integer.")
    end
    
    # Check for `n`
    if typeof(yaml_config["n"]) != Int || yaml_config["n"] < 0
        error("Invalid value for 'n': must be a non-negative integer.")
    end

    # Check for `n_points`
    if yaml_config["n_points"] != 0 && (typeof(yaml_config["n_points"]) != Vector{Int} || length(yaml_config["n_points"]) != yaml_config["k"])
        error("Invalid configuration for 'n_points': must be 0 or a list of $(yaml_config["k"]) integers, instead found $(yaml_config["n_points"])")
    end

    # Check for `ub_n_points` (optional, used only when `n_points == 0`)
    if haskey(yaml_config, "ub_n_points")
        if typeof(yaml_config["ub_n_points"]) != Int || yaml_config["ub_n_points"] < 1
            error("Invalid value for 'ub_n_points': must be a positive integer.")
        end
    end

    # Check for `vertex_sampling` (optional)
    vertex_sampling = get(yaml_config, "vertex_sampling", "box_uniform")
    if !(vertex_sampling isa String) || !(vertex_sampling in ["box_uniform", "sphere", "ellipsoid"])
        error("Invalid value for 'vertex_sampling': must be 'box_uniform', 'sphere', or 'ellipsoid'.")
    end

    # Check for `sphere_radius_factor` (optional, used only when `vertex_sampling == \"sphere\"` or \"ellipsoid\"`)
    if vertex_sampling in ["sphere", "ellipsoid"] || haskey(yaml_config, "sphere_radius_factor")
        sphere_radius_factor = get(yaml_config, "sphere_radius_factor", 0.95)
        if !(sphere_radius_factor isa Float64) || sphere_radius_factor <= 0.0 || sphere_radius_factor >= 1.0
            error("Invalid value for 'sphere_radius_factor': must be a Float64 strictly between 0 and 1.")
        end
    end

    # Check for `target_tolerance`
    if typeof(yaml_config["target_tolerance"]) != Float64 || yaml_config["target_tolerance"] < 0
        error("Invalid value for 'target_tolerance': must be a non-negative float.")
    end

    # Check for `target_tolerance_opt`
    if typeof(yaml_config["target_tolerance_opt"]) != Float64 || yaml_config["target_tolerance_opt"] < 0
        error("Invalid value for 'target_tolerance_opt': must be a non-negative float.")
    end

    # Check for `max_iterations`
    if typeof(yaml_config["max_iterations"]) != Int || yaml_config["max_iterations"] < 1
        error("Invalid value for 'max_iterations': must be a positive integer.")
    end

    # Check for `max_iterations_opt`
    if typeof(yaml_config["max_iterations_opt"]) != Int || yaml_config["max_iterations_opt"] < 1
        error("Invalid value for 'max_iterations_opt': must be a positive integer.")
    end

    # Check for `max_print_iterations`
    if typeof(yaml_config["max_print_iterations"]) != Int || yaml_config["max_print_iterations"] < 1
        error("Invalid value for 'max_print_iterations': must be a positive integer.")
    end

    # Check for `seed`
    if typeof(yaml_config["seed"]) != Int || yaml_config["seed"] < 0
        error("Invalid value for 'seed': must be a positive integer.")
    end

    # Check for `cvxhflag`
    if typeof(yaml_config["cvxhflag"]) != Bool
        error("Invalid value for 'cvxhflag': must be a boolean.")
    end
    # Check for `convex_hull_backend` (optional; only used when `cvxhflag == true`)
    if haskey(yaml_config, "convex_hull_backend")
        backend = yaml_config["convex_hull_backend"]
        if !(backend isa String) || !(backend in ("matrix", "vector"))
            error("Invalid value for 'convex_hull_backend': must be either 'matrix' or 'vector'.")
        end
    end
    # Check for `verbose` (optional)
    if haskey(yaml_config, "verbose")
        if !(yaml_config["verbose"] isa Bool)
            error("Invalid value for 'verbose': must be a boolean.")
        end
    end
    # Check for `matrix_lmo_cache_cap` (optional)
    if haskey(yaml_config, "matrix_lmo_cache_cap")
        cap = yaml_config["matrix_lmo_cache_cap"]
        if !(cap isa Int) || cap < -1
            error("Invalid value for 'matrix_lmo_cache_cap': must be an integer >= -1.")
        end
    end
    # Check for `matrix_lmo_use_optimized_search` (optional)
    if haskey(yaml_config, "matrix_lmo_use_optimized_search")
        if !(yaml_config["matrix_lmo_use_optimized_search"] isa Bool)
            error("Invalid value for 'matrix_lmo_use_optimized_search': must be a boolean.")
        end
    end

    # Check for `intersection`
    if !haskey(yaml_config, "intersection") || !(yaml_config["intersection"] isa Dict)
        error("Missing or invalid 'intersection' section: must be a mapping with 'anchor' and 'reference_point'.")
    end
    intersection_cfg = yaml_config["intersection"]
    # Check for `intersection.anchor`
    if !haskey(intersection_cfg, "anchor") || !(intersection_cfg["anchor"] isa String)
        error("Invalid or missing 'intersection.anchor': must be a string.")
    end
    if !(intersection_cfg["anchor"] in ["p1_center", "p1_vertex", "p1_random", "origin", "global_mean", "random", "p1_p2_segment"])
        error("Invalid value for 'intersection.anchor': must be one of 'p1_center', 'p1_vertex', 'p1_random', 'origin', 'global_mean', 'random', 'p1_p2_segment'.")
    end
    # Check for `intersection.anchor_t` (optional, used only for some anchors)
    if haskey(intersection_cfg, "anchor_t")
        if !(intersection_cfg["anchor_t"] isa Float64) || intersection_cfg["anchor_t"] < 0.0
            error("Invalid value for 'intersection.anchor_t': must be a non-negative float.")
        end
    end
    # Check for `intersection.reference_point`
    if !haskey(intersection_cfg, "reference_point") || !(intersection_cfg["reference_point"] isa String)
        error("Invalid or missing 'intersection.reference_point': must be a string ('center' or 'vertex').")
    end
    if !(intersection_cfg["reference_point"] in ["center", "vertex"])
        error("Invalid value for 'intersection.reference_point': must be 'center' or 'vertex'.")
    end

    # Check for `stepsize_strategy`
    if typeof(yaml_config["stepsize_strategy"]) != Int || !(yaml_config["stepsize_strategy"] in (0, 1))
        error("Invalid value for 'stepsize_strategy': must be an integer in {0, 1}.")
    end
end

# Function to print the config parameters
function print_config(config::Config)
    
    println()
    println("Configuration Parameters:")
    println("  Number of polytopes (k): ", config.k)
    println("  Polytope dimension (n): ", config.n)
    println("  Number of points (n_points): ", config.n_points)
    println("  n_points upper bound multiplier (ub_n_points): ", config.ub_n_points)
    println("  Vertex sampling strategy (vertex_sampling): ", config.vertex_sampling)
    println("  Sphere/ellipsoid factor (sphere_radius_factor): ", config.sphere_radius_factor)
    println("  Epsilon-optimality threshold (target_tolerance): ", config.target_tolerance)
    println("  Epsilon-optimality threshold to compute optimal sols (target_tolerance_opt): ", config.target_tolerance_opt)
    println("  Number of FW iterations (max_iterations): ", config.max_iterations)
    println("  Number of FW iterations to compute optimal solutions (max_iterations_opt): ", config.max_iterations_opt)
    println("  How often FW iteration log is printed to screen (max_print_iterations): ", config.max_print_iterations)
    println("  Seed for reproducibility (seed): ", config.seed)
    println("  Verbose diagnostics (verbose): ", config.verbose)
    println("  Use convex-hull LMOs or MathOptLMO (cvxhflag): ", config.cvxhflag)
    println("  Convex-hull backend (convex_hull_backend): ", config.convex_hull_backend)
    println("  MatrixConvexHullLMO cache cap (matrix_lmo_cache_cap): ", config.matrix_lmo_cache_cap)
    println("  MatrixConvexHullLMO optimized search (matrix_lmo_use_optimized_search): ", config.matrix_lmo_use_optimized_search)
    println("  Intersection anchor: ", config.intersection_anchor)
    println("  Intersection reference point: ", config.intersection_reference_point)
    println("  Intersection anchor_t: ", config.intersection_anchor_t)
    println("  Stepsize strategy (`0`: SafeGoldenratio (goldensearch), `1`: Shortstep(L=1)): ", config.stepsize_strategy)
    println()
end

function get_stepsize_strategy(stepsize_strategy::Int, L::T) where T
    
    if stepsize_strategy == 0
        return SafeGoldenratio(1e-09)           # same line search, safer gamma reconstruction
    elseif stepsize_strategy == 1
        return FrankWolfe.Shortstep(L)          # short step with given L
    else
        error("Invalid stepsize strategy type")
    end
end

# turns Config objectinto dictionary that YAML.jl can serialise
function todict(config::Config)
    # Note: we explicitly build the dictionary so that the YAML layout matches examples/config.yml,
    # in particular nesting intersection settings under the "intersection" key and not exposing anc_flag.
    return Dict{String, Any}(
        "k" => config.k,
        "n" => config.n,
        "n_points" => config.n_points,
        "ub_n_points" => config.ub_n_points,
        "vertex_sampling" => config.vertex_sampling,
        "sphere_radius_factor" => config.sphere_radius_factor,
        "target_tolerance" => config.target_tolerance,
        "target_tolerance_opt" => config.target_tolerance_opt,
        "max_iterations" => config.max_iterations,
        "max_iterations_opt" => config.max_iterations_opt,
        "max_print_iterations" => config.max_print_iterations,
        "seed" => config.seed,
        "verbose" => config.verbose,
        "cvxhflag" => config.cvxhflag,
        "convex_hull_backend" => config.convex_hull_backend,
        "matrix_lmo_cache_cap" => config.matrix_lmo_cache_cap,
        "matrix_lmo_use_optimized_search" => config.matrix_lmo_use_optimized_search,
        "intersection" => Dict(
            "anchor" => config.intersection_anchor,
            "reference_point" => config.intersection_reference_point,
            "anchor_t" => config.intersection_anchor_t,
        ),
        "stepsize_strategy" => config.stepsize_strategy,
    )
end
    
# Serialize `config` to YAML and store it at `config_filename`
function write_config(config::Config, config_filename::AbstractString)
    
    # make sure the destination directory exists
    mkpath(dirname(config_filename))
    YAML.write_file(config_filename, todict(config))
    return abspath(config_filename)
end

function generate_n_points(n::Integer, k::Integer, seed::Integer; ub_n_points::Integer=10)
    # Generate a list of k random integers in {n+1, ..., n*ub_n_points}
    Random.seed!(seed)  # Set the seed for reproducibility
    upper_bound = max(n + 1, n * ub_n_points)
    return [rand((n + 1):upper_bound) for _ in 1:k]
end
