# `config.jl`
struct Config
    # Number of polytopes
    k::Int
    # Dimension of the subspace in which each polytope lies
    n::Int
    # list of `k` integers, each indicating the n. of points to be used to generate intersecting polytopes 
    n_points::Vector{Int}
    # epsilon-optimality threshold
    target_tolerance::Float64
    # Number of FW iterations
    max_iterations::Int
    # Set the seed for reproducibility
    seed::Int
end

# Constructor with default values
function Config()
    return Config(2, 2, [5, 5], 1e-6, 10000, 42)
end

# Constructor from YAML file
function Config(yaml_file::String)
    # Load and validate configuration from YAML file
    config = YAML.load(open(yaml_file))
    
    # Validation checks (for YAML file dict)
    validate_config(config)

    return Config(config["k"], config["n"], config["n_points"], config["target_tolerance"], config["max_iterations"], config["seed"])
end

# Validation checks (for YAML file dict)
function validate_config(config)
    # Check for `k`
    if typeof(config["k"]) != Int || config["k"] < 1
        error("Invalid value for 'k': must be a positive integer.")
    end
    
    # Check for `n`
    if typeof(config["n"]) != Int || config["n"] < 0
        error("Invalid value for 'n': must be a non-negative integer.")
    end
    
    # Assuming 'config' is a dictionary and 'k' is a key in the dictionary
    if typeof(config["n_points"]) != Vector{Int} || length(config["n_points"]) != config["k"]
        error("Invalid configuration for 'n_points': must be a list of $(config["k"]) integers.")
    end

    # Check for `target_tolerance`
    if typeof(config["target_tolerance"]) != Float64 || config["target_tolerance"] < 0
        error("Invalid value for 'target_tolerance': must be a non-negative float.")
    end

    # Check for `max_iterations`
    if typeof(config["max_iterations"]) != Int || config["max_iterations"] < 1
        error("Invalid value for 'max_iterations': must be a positive integer.")
    end

    # Check for `seed`
    if typeof(config["seed"]) != Int || config["seed"] < 0
        error("Invalid value for 'seed': must be a positive integer.")
    end
end
