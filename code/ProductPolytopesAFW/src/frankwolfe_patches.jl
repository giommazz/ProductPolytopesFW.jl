# `frankwolfe_patches.jl`
# Small opt-in line-search patches for FrankWolfe.jl behavior.


# -----------------------------------------------------------------------------
# Line Search Extensions
# -----------------------------------------------------------------------------
"""
    SafeGoldenratio(tol=1e-7)

Like `FrankWolfe.Goldenratio` but with numerically safer reconstruction of final step size `γ`.

Why:
- `Goldenratio` computes bracket `[left, right]` in primal space, then reconstructs `γ` via a
  single-coordinate division using *1st* nonzero entry of direction `d`.
  When `d[i]` is tiny (structured atoms, sparse directions, etc.), division can be
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

    # Restrict segment of search to [x, y] where y = x - gamma_max*d
    @. workspace.y = x - gamma_max * d
    @. workspace.left = x
    @. workspace.right = workspace.y

    # Compute directional derivatives along search direction at endpoints.
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
        # Maintain a bracket [left, right], replace worse side based on two interior probes (classic GoldenRatio)
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

    # Reconstruct gamma from final bracket using projection formula:
    #   x_min = (left + right)/2, and we want gamma s.t. x_min ≈ x - gamma*d.
    # More stable: error in GoldenRatio scales as 1/|dᵢ|, error in SafeGoldenRatio scales as 1/||d||₂, and
    #   |dᵢ| = |<d, eᵢ>| ≤ ||d||₂⋅||eᵢ||₂ = ||d||₂
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
    # Final guard rail: clamp. Rounding can produce gamma that is tiny negative or
    # slightly above `gamma_max` even if true minimizer is ∈ [0, gamma_max].
    return clamp(gamma, zero(gamma), gamma_max)
end

Base.print(io::IO, ::SafeGoldenratio) = print(io, "SafeGoldenratio")