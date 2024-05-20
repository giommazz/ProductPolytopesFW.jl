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
    # Number of FW iterations
    max_iterations::Int
    # Set the seed for reproducibility
    seed::Int
end

# Constructor with default values
function Config()
    c = Config(2, 2, [5, 5], 1e-6, 10000, 42)
    print_config(c)
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
        # Generate a list of random integers in [3, 200] of length config["k"]
        Random.seed!(config["seed"])  # Set the seed for reproducibility
        n_points = [rand(config["n"]:config["n"]*2 + 1) for _ in 1:config["k"]]
    end
    c = Config(config["k"], config["n"], n_points, config["target_tolerance"], config["max_iterations"], config["seed"])
    print_config(c)
    return c 
end

# Function to update Config with any number of arguments
function update_config(config::Config; kwargs...)
    # `get(kwargs, :key, default_value)` retrieves value for `:key` from `kwargs` if it exists, o/w returns default_value
    k = get(kwargs, :k, config.k)
    n = get(kwargs, :n, config.n)
    n_points = get(kwargs, :n_points, config.n_points)
    target_tolerance = get(kwargs, :target_tolerance, config.target_tolerance)
    max_iterations = get(kwargs, :max_iterations, config.max_iterations)
    seed = get(kwargs, :seed, config.seed)

    # Return a new Config object with updated values
    c = Config(k, n, n_points, target_tolerance, max_iterations, seed)
    print_config(c)
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

    # Check for `max_iterations`
    if typeof(yaml_config["max_iterations"]) != Int || yaml_config["max_iterations"] < 1
        error("Invalid value for 'max_iterations': must be a positive integer.")
    end

    # Check for `seed`
    if typeof(yaml_config["seed"]) != Int || yaml_config["seed"] < 0
        error("Invalid value for 'seed': must be a positive integer.")
    end
end

# Function to print the config parameters
function print_config(config::Config)
    println("Configuration Parameters:")
    println("  Number of polytopes (k): ", config.k)
    println("  Dimension of the subspace (n): ", config.n)
    println("  Number of points (n_points): ", config.n_points)
    println("  Epsilon-optimality threshold (target_tolerance): ", config.target_tolerance)
    println("  Number of FW iterations (max_iterations): ", config.max_iterations)
    println("  Seed for reproducibility (seed): ", config.seed)
end