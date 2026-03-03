# `frankwolfe_patches.jl`
# Small opt-in line-search patches for FrankWolfe.jl behavior.


# -----------------------------------------------------------------------------
# Line Search Extensions
# -----------------------------------------------------------------------------
"""
    SafeGoldenratio(tol=1e-7)

Golden-section line search like `FrankWolfe.Goldenratio`, but with a numerically safer
reconstruction of the final step size `γ`.

Why:
- In FrankWolfe.jl v0.6.3, `Goldenratio` computes bracket `[left, right]` in primal space, then 
  reconstructs `γ` via a single-coordinate division using *1st* nonzero entry of direction `d`.
- When `d[i]` is tiny (common with structured atoms / sparse directions), division can be
  ill-conditioned, producing inaccurate `γ` (negative in worst case) and nonmonotone primal gap.

`SafeGoldenratio` computes `γ` using all coordinates via a projection formula and clamps it to
`[0, gamma_max]`.
"""
struct SafeGoldenratio{T} <: FrankWolfe.LineSearchMethod
    tol::T
end

SafeGoldenratio() = SafeGoldenratio(1e-7)

struct SafeGoldenratioWorkspace{XT,GT}
    y::XT
    left::XT
    right::XT
    new_vec::XT
    probe::XT
    gradient::GT
end

function FrankWolfe.build_linesearch_workspace(::SafeGoldenratio, x::XT, gradient::GT) where {XT,GT}
    return SafeGoldenratioWorkspace{XT,GT}(
        similar(x),
        similar(x),
        similar(x),
        similar(x),
        similar(x),
        similar(gradient),
    )
end

function FrankWolfe.perform_line_search(
    line_search::SafeGoldenratio,
    _,
    f,
    grad!,
    gradient,
    x,
    d,
    gamma_max,
    workspace::SafeGoldenratioWorkspace,
    memory_mode,
)
    # Intentionally mirrors FrankWolfe.jl's `Goldenratio` bracketing logic.
    # The "patch" is only how we reconstruct the final scalar stepsize `gamma` at the end.
    #
    # Why we care: FW direction `d` can have some coordinates that are *very* small (due to
    # cancellations/symmetries) and reconstructing `gamma` via a single coordinate division
    #   (x[i]-x_min[i]) / d[i]
    # can become numerically unstable if chosen d[i] is tiny.

    # Restrict segment of search to [x, y] where y = x - gamma_max*d
    @. workspace.y = x - gamma_max * d
    @. workspace.left = x
    @. workspace.right = workspace.y

    # Directional derivatives along search direction at endpoints.
    # If they do not change sign, minimizer over segment is at an endpoint.
    dgx = dot(d, gradient)
    grad!(workspace.gradient, workspace.y)
    dgy = dot(d, workspace.gradient)
    # If minimum is at an endpoint
    if dgx * dgy >= 0
        return f(workspace.y) <= f(x) ? gamma_max : zero(eltype(d))
    end

    # Apply golden-section method to segment
    gold = (1 + sqrt(5)) / 2
    improv = Inf
    while improv > line_search.tol
        # Classic golden-section shrink: maintain a bracket [left, right], replace "worse"
        # side based on two interior probes.
        f_old_left = f(workspace.left)
        f_old_right = f(workspace.right)
        @. workspace.new_vec = workspace.left + (workspace.right - workspace.left) / (1 + gold)
        @. workspace.probe = workspace.new_vec + (workspace.right - workspace.new_vec) / 2

        if f(workspace.probe) <= f(workspace.new_vec)
            workspace.left .= workspace.new_vec
        else
            workspace.right .= workspace.probe
        end
        # Compute stopping criterion quantity
        improv = norm(f(workspace.right) - f_old_right) + norm(f(workspace.left) - f_old_left)
    end

    # Reconstruct gamma from final bracket using stable formula:
    #   x_min = (left + right)/2, and we want gamma s.t. x_min ≈ x - gamma*d.
    #
    # Arithmetically, `x_min` lies on the same line segment, so *any* nonzero coordinate
    # would give same gamma via (x[i]-x_min[i]) / d[i]. Numerically, though, picking a single
    # coordinate can be bad when that d[i] is tiny → projection formula below uses all 
    # coordinates and is better conditioned
    den = dot(d, d)
    if den <= eps(float(den))
        return zero(eltype(d))
    end
    num = zero(den)
    @inbounds for i in eachindex(d)
        x_min_i = (workspace.left[i] + workspace.right[i]) / 2
        num += (x[i] - x_min_i) * d[i]
    end
    gamma = num / den
    if !isfinite(gamma)
        return zero(eltype(d))
    end
    # Final guard rail: clamp. Rounding can produce a tiny negative gamma or a gamma
    # slightly above `gamma_max` even if true minimizer is ∈ [0, gamma_max].
    return clamp(gamma, zero(gamma), gamma_max)
end

Base.print(io::IO, ::SafeGoldenratio) = print(io, "SafeGoldenratio")