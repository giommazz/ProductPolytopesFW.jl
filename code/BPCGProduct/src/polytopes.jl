# `polytopes.jl`
const MOI = MathOptInterface

# (Multiple Dispatch) Generate a random convex polytope A in ℝⁿ with m constraints
function generate_polytope(m, config)
    Random.seed!(config.seed)
    model_A = Model()
    @variable(model_A, x[1:config.n] ≥ 0)  # Variable representing the polytope

    # Generate random integer coefficients for matrix A and vector b
    A = [rand(-10:10) for _ in 1:m, _ in 1:config.n]
    b = [rand(0:10) for _ in 1:m]  # Example: vector with random positive integers between 0 and 10

    # Add constraints to the model
    for i in 1:m
        @constraint(model_A, dot(A[i, :], x) ≤ b[i])
    end
    return model_A, x
end
# (Multiple Dispatch) Generate a pyramid-like polytope B in ℝⁿ    
function generate_polytope(config)
    Random.seed!(config.seed)
    model_B = Model()
    @variable(model_B, x[1:config.n])  # Variable representing the polytope

    # Tip of the pyramid (random integer coefficient)
    @constraint(model_B, x[1] ≥ rand(0:10))

    # Base constraints for the pyramid
    for i in 2:config.n
        @constraint(model_B, x[i] ≥ -rand(1:10) * x[1])
        @constraint(model_B, x[i] ≤ rand(1:10) * x[1])
    end
    return model_B, x
end

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

# # Function to plot a 2D polytope defined by constraints in a JuMP model
# function plot_polytope_2d(model::Model, x_vars)
#     # Retrieve constraints with AffExpr functions and LessThan set types
#     less_than_constraints = all_constraints(model, AffExpr, MOI.LessThan{Float64})

#     # Initialize a plot with axis limits for visualization
#     plot(title="2D Polytope", xlim=(-2, 2), ylim=(-2, 2), legend=false)

#     # Iterate through the retrieved constraints and plot them
#     for constraint in less_than_constraints
#         expr = constraint.func
#         d = constraint.set.upper

#         # Extract coefficients from the affine expression (dot product)
#         c = JuMP.coefficients(expr, x_vars)

#         # Check if the equation has non-zero coefficients for plotting
#         if length(c) == 2 && c[2] != 0
#             # Rearrange to y = mx + c
#             m = -c[1] / c[2]
#             intercept = d / c[2]
#             plot!(x -> m * x + intercept, label=false)
#         end
#     end

#     display(plot)
# end