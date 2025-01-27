# `config.jl`
struct Config
    
    # Number of polytopes
    k::Int
    # Dimension of the subspace in which each polytope lies
    n::Int
    # Either 0 or a list of `k` integers, each indicating the number of points to be used to generate intersecting polytopes
    n_points::Union{Int, Vector{Int}}
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
    # Use FrankWolfe.ConvexHullOracle LMOs (true) or FrankWolfe.MathOptLMO LMOs (false)
    cvxhflag::Bool
    # Intersect polytopes close to the first polytope's analytic center (true) or a vertex
    anc_flag::Bool
    # stepsize strategy: `0` is line search; `1` uses short-step with given L=2 (specified in `product_algorithms.jl`)
    stepsize_strategy::Int 
end

# Constructor with default values
function Config()
    
    c = Config(2, 2, [5, 5], 1e-6, 1e-08, 1000, 5000, 100, 42, true, false, 0)
    
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
    if n_points == 0
        # Generate a list of random integers in [n, 2n+1] of length config["k"]
        Random.seed!(config["seed"])  # Set the seed for reproducibility
        n_points = [rand(config["n"]:config["n"]*2 + 1) for _ in 1:config["k"]]
    end
    c = Config(
            config["k"],
            config["n"],
            n_points,
            config["target_tolerance"],
            config["target_tolerance_opt"],
            config["max_iterations"],
            config["max_iterations_opt"],
            config["max_print_iterations"],
            config["seed"],
            config["cvxhflag"],
            config["anc_flag"],
            config["stepsize_strategy"]
            )
    return c 
end

# Function to update Config with any number of arguments
function modify_config(config::Config; kwargs...)
    
    # `get(kwargs, :key, default_value)` retrieves value for `:key` from `kwargs` if it exists, o/w returns default_value
    k = get(kwargs, :k, config.k)
    n = get(kwargs, :n, config.n)
    n_points = get(kwargs, :n_points, config.n_points)
    target_tolerance = get(kwargs, :target_tolerance, config.target_tolerance)
    target_tolerance_opt = get(kwargs, :target_tolerance_opt, config.target_tolerance_opt)
    max_iterations = get(kwargs, :max_iterations, config.max_iterations)
    max_iterations_opt = get(kwargs, :max_iterations_opt, config.max_iterations_opt)
    max_print_iterations = get(kwargs, :max_iterations, config.max_print_iterations)
    seed = get(kwargs, :seed, config.seed)
    cvxhflag = get(kwargs, :cvxhflag, config.cvxhflag)
    anc_flag = get(kwargs, :anc_flag, config.anc_flag)
    stepsize_strategy = get(kwargs, :stepsize_strategy, config.stepsize_strategy)

    # Return a new Config object with updated values
    c = Config(k, n, n_points, target_tolerance, target_tolerance_opt, max_iterations, max_iterations_opt, max_print_iterations, seed, cvxhflag, anc_flag, stepsize_strategy)
    
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
        error("Invalid configuration for 'n_points': must be 0 or a list of $(yaml_config["k"]) integers.")
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

    # Check for `anc_flag`
    if typeof(yaml_config["anc_flag"]) != Bool
        error("Invalid value for 'anc_flag': must be a boolean.")
    end

    # Check for `stepsize_strategy`
    if typeof(yaml_config["stepsize_strategy"]) != Int || yaml_config["stepsize_strategy"] < 0 || yaml_config["stepsize_strategy"] > 1
        error("Invalid value for 'stepsize_strategy': must be a positive integer in {0, 1}.")
    end
end

# Function to print the config parameters
function print_config(config::Config)
    
    println()
    println("Configuration Parameters:")
    println("  Number of polytopes (k): ", config.k)
    println("  Dimension of the subspace (n): ", config.n)
    println("  Number of points (n_points): ", config.n_points)
    println("  Epsilon-optimality threshold (target_tolerance): ", config.target_tolerance)
    println("  Epsilon-optimality threshold to compute optimal sols (target_tolerance_opt): ", config.target_tolerance_opt)
    println("  Number of FW iterations (max_iterations): ", config.max_iterations)
    println("  Number of FW iterations to compute optimal solutions (max_iterations_opt): ", config.max_iterations_opt)
    println("  How often FW iteration log is printed to screen (max_print_iterations): ", config.max_print_iterations)
    println("  Seed for reproducibility (seed): ", config.seed)
    println("  Use FW's ConvexHullOracle LMOs or MathOptLMO (cvxhflag): ", config.cvxhflag)
    println("  Intersect polytopes close to first polytope's analytic center (true) or vertex (anc_flag): ", config.anc_flag)
    println("  Stepsize strategy (`0` is line search; `1` uses short-step with L=2 specified in `product_algorithms.jl`): ", config.stepsize_strategy)
    println()
end

function get_stepsize_strategy(stepsize_strategy::Int, L::T) where T
    if stepsize_strategy == 0
        return FrankWolfe.Goldenratio()         # simple line search
    elseif stepsize_strategy == 1
        return FrankWolfe.Shortstep(L)          # short step with given L
    else
        error("Invalid stepsize strategy type")
    end
end