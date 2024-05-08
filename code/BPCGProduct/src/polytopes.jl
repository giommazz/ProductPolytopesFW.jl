# `polytopes.jl`
const MOI = MathOptInterface

# GENERATING RANDOM INTEGER NUMBERS
# (Multiple Dispatch) Generate a random convex polytope A in ℝⁿ with m constraints
function generate_polytope(m, config)
    Random.seed!(config.seed)
    model_A = Model(GLPK.Optimizer)
    @variable(model_A, x[1:config.n] ≥ 0)  # Variable representing the polytope

    # Generate random integer coefficients for matrix A and vector b
    A = [rand(-10:10) for _ in 1:m, _ in 1:config.n]
    b = [rand(0:10) for _ in 1:m]  # Example: vector with random positive integers between 0 and 10

    # Add constraints to the model
    for i in 1:m
        @constraint(model_A, dot(A[i, :], x) ≤ b[i])
    end
    return model_A, x, A, b
end
# Solve for a vertex of polytope A by maximizing an objective function
function find_vertex_in_polytope(model_A, x, direction)
    @objective(model_A, Max, dot(direction, x))
    optimize!(model_A)
    return value.(x)
end
# Generate a simplex-like pointy polytope in R^n
function generate_simplex(config)
    model = Model()
    @variable(model, x[1:config.n])  # Variable representing the polytope

    # Define n+1 constraints representing a regular simplex
    for i in 1:config.n
        @constraint(model, x[i] >= 0)
    end

    @constraint(model, sum(x) <= 1)

    return model, x
end


function setup_translated_polytope_B(config, vertex_A)
    model_B, x_B = generate_simplex(config)

    # Apply translation to polytope B (shift constraints)
    translated_constraints = @variable(model_B, translated_x[1:config.n])

    for i in 1:config.n
        @constraint(model_B, translated_x[i] == x_B[i] + vertex_A[i])
    end

    return model_B, translated_constraints
end





# Generate a shifted simplex-like polytope in R^n
function generate_shifted_simplex_polytope(config, offset)
    model = Model()
    @variable(model, x[1:config.n])  # Variables representing the polytope

    # 1. Create the shifted variables explicitly using an indexed for loop and applying an offset to each
    # x_shifted = [@variable(model, x[i] + offset[i]) for i in 1:config.n]
    # 2. Apply an offset to all variables
    # @expression(model, x_shifted[i in 1:config.n] == x[i] + offset[i])

    # 3.1. Create new variables shifted by the provided offset
    x_shifted = [@variable(model) for _ in 1:config.n]
    
    println("°°°°°°°°°°°°°°°°°°°°°°°", offset)

    # 3.2. Apply an offset to all shifted variables through constraints
    for i in 1:config.n
        @constraint(model, x_shifted[i] == x[i] + offset[i])
    end

    # Define n+1 constraints representing a regular simplex
    for i in 1:config.n
        @constraint(model, x_shifted[i] >= 0)
    end

    @constraint(model, sum(x_shifted) <= 1)
    println("°°°°°°°°°°°°°°°°°°°°°°°", x)
    println("°°°°°°°°°°°°°°°°°°°°°°°", x_shifted)
    return model, x
end




# GENERATING RANDOM FLOAT NUMBERS
# # (Multiple Dispatch) Generate a random convex polytope A in ℝⁿ with m constraints
# function generate_polytope(m, config)
#     Random.seed!(config.seed)
#     model_A = Model()
#     @variable(model_A, x[1:config.n] ≥ 0)  # Variable representing the polytope

#     # Randomly generate coefficient matrix and vector
#     A = rand(m, config.n)
#     b = rand(m)

#     # Add constraints to the model
#     for i in 1:m
#         @constraint(model_A, dot(A[i, :], x) ≤ b[i])
#     end
#     return model_A, x
# end
# # (Multiple Dispatch) Generate a pyramid-like polytope B in ℝⁿ    
# function generate_polytope(config)
#     model_B = Model()
#     @variable(model_B, x[1:config.n])  # Variable representing the polytope

#     # Tip of the pyramid
#     @constraint(model_B, x[1] >= 0)

#     # Base constraints for the pyramid
#     for i in 2:config.n
#         @constraint(model_B, x[i] >= -x[1])
#         @constraint(model_B, x[i] <= x[1])
#     end
#     return model_B, x
# end

# Function to generate a polytope with a given JuMP model
function polytope_from_jump(model)
    # Convert JuMP model to a Polyhedra polytope
    polytope = polyhedron(model, CDDLib.Library(:exact))
    polymesh = Polyhedra.Mesh(polytope)
    fig = Makie.mesh(polymesh, color=:blue)
    display(fig)
    return polytope, fig
end


TODO: 
- decide n-dimensional box A
- decide another n-dimensional box B, such that the two boxes don't intersect
- generate m vertices in A, and compute convex hull with Polyhedra.jl --> obtain polytope X
- do the same with B --> obtain polytope Y
- transform X into JuMP
- optimize in random direction d over X -> x1
- optimize in direction -d over Y -> y1
- now take each vertex of Y and move it in the direction x1 - y1. so for each vertex y in Y: do y += x1 - y1
- now the two polytopes touch in at least one point and the intersection is not empty
