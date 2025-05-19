# `compute_intersection_custom_instances.jl`
using ProductPolytopesAFW
using FrankWolfe
using LinearAlgebra

# EXAMPLE WITH HYPERCUBE
struct HypercubeLMO <: FrankWolfe.LinearMinimizationOracle
    n::Int
end

function FrankWolfe.compute_extreme_point(lmo::HypercubeLMO, direction::AbstractVector; v=nothing)
    n = lmo.n
    result = similar(direction)
    for i in 1:length(direction)
        result[i] = direction[i] < 0 ? -1.0 : 1.0/n
    end
    return result
end

# Objective function
function f(x)
    return sum(abs2, x)  # Example objective function
end

# Gradient of the objective function
function grad!(storage, x)
    storage .= 2x  # In-place gradient computation
end

n = 10

# Initialize LMO
lmo = HypercubeLMO(n)

# Initial point
x0 = randn(n * n)

# Run Frank-Wolfe algorithm
x_opt, _ = frank_wolfe(f, grad!, lmo, x0; max_iteration=100, verbose=true)


# EXAMPLE WITH SPECTRAPLEX

# Define the objective function
function f(X)
    return sum(X.^2)  # Example objective function
end

# Gradient of the objective function
function grad!(storage, X)
    storage .= 2X  # In-place gradient computation
end

# Initialize Spectraplex LMO
n = 10  # example dimension
radius = 1.0
ensure_symmetry = true
gradient_container = zeros(n, n)
lmo = FrankWolfe.SpectraplexLMO(radius, gradient_container, ensure_symmetry, n)

# Initial point
X0 = randn(n, n)
X0 = X0 * X0'  # To ensure it's PSD
X0 ./= LinearAlgebra.tr(X0)  # Normalize to make trace 1

# Run Frank-Wolfe algorithm
X_opt, _ = frank_wolfe(f, grad!, lmo, X0; max_iteration=100, verbose=true)