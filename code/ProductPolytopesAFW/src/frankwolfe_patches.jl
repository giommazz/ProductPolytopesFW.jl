# `frankwolfe_patches.jl`
#
# Small opt-in monkey patches / overrides for FrankWolfe.jl behavior.
#
# These are meant for controlled experiments and debugging. Prefer keeping them disabled by default.

const FW_WEIGHT_PURGE_DEFAULT_OVERRIDE = Ref{Union{Nothing,Float64}}(nothing)

"""
    set_fw_weight_purge_default_override!(value::Union{Nothing,Float64}) -> Union{Nothing,Float64}

Override the value returned by `FrankWolfe.weight_purge_threshold_default(Float64)`.

Why this exists:
- In FrankWolfe.jl, internal active-set cleanups during AFW iterations call
  `active_set_update!(...; weight_purge_threshold=weight_purge_threshold_default(R))`.
  Since `away_frank_wolfe` does not pass your `weight_purge_threshold` down to those calls,
  changing `weight_purge_threshold` mostly affects the *final* cleanup, not the in-run ones.
- Overriding the default for `Float64` lets you test whether those internal purges are the culprit.

Pass `nothing` to restore FrankWolfe's original default (≈ 1e-12 for Float64).
"""
function set_fw_weight_purge_default_override!(value::Union{Nothing,Float64})
    if value !== nothing && value < 0.0
        error("`value` must be ≥ 0 or `nothing`, got $value")
    end
    FW_WEIGHT_PURGE_DEFAULT_OVERRIDE[] = value
    return value
end

"""
    fw_weight_purge_default_override() -> Union{Nothing,Float64}

Return the currently configured override for `FrankWolfe.weight_purge_threshold_default(Float64)`.
"""
fw_weight_purge_default_override() = FW_WEIGHT_PURGE_DEFAULT_OVERRIDE[]

_fw_weight_purge_default_float64() = sqrt(eps(Float64) * Base.rtoldefault(Float64))

function FrankWolfe.weight_purge_threshold_default(::Type{Float64})
    v = FW_WEIGHT_PURGE_DEFAULT_OVERRIDE[]
    return v === nothing ? _fw_weight_purge_default_float64() : v
end


# -----------------------------------------------------------------------------
# Line Search Extensions
# -----------------------------------------------------------------------------

"""
    SafeGoldenratio(tol=1e-7)

Golden-section line search like `FrankWolfe.Goldenratio`, but with a numerically safer
reconstruction of the final step size `γ`.

Why this exists:
- In FrankWolfe.jl v0.6.3, `Goldenratio` computes a bracket `[left, right]` in *primal space*
  and then reconstructs `γ` via a single coordinate division using the *first* nonzero entry
  of the direction `d`.
- When `d[i]` is tiny (common with structured atoms / sparse directions), that division can be
  ill-conditioned, producing an inaccurate `γ` (and in the worst case a negative one).

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

    dgx = dot(d, gradient)
    grad!(workspace.gradient, workspace.y)
    dgy = dot(d, workspace.gradient)

    # If the minimum is at an endpoint
    if dgx * dgy >= 0
        return f(workspace.y) <= f(x) ? gamma_max : zero(eltype(d))
    end

    # Apply golden-section method to the segment
    gold = (1 + sqrt(5)) / 2
    improv = Inf
    while improv > line_search.tol
        f_old_left = f(workspace.left)
        f_old_right = f(workspace.right)

        @. workspace.new_vec = workspace.left + (workspace.right - workspace.left) / (1 + gold)
        @. workspace.probe = workspace.new_vec + (workspace.right - workspace.new_vec) / 2

        if f(workspace.probe) <= f(workspace.new_vec)
            workspace.left .= workspace.new_vec
        else
            workspace.right .= workspace.probe
        end

        improv = norm(f(workspace.right) - f_old_right) + norm(f(workspace.left) - f_old_left)
    end

    # Reconstruct gamma from the final bracket using a stable formula:
    # x_min = (left + right)/2, and we want gamma s.t. x_min ≈ x - gamma*d.
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
    return clamp(gamma, zero(gamma), gamma_max)
end

Base.print(io::IO, ::SafeGoldenratio) = print(io, "SafeGoldenratio")


"""
    QuadraticExactLineSearch()

Exact line search for *homogeneous quadratic* objectives of the form `f(z) = 0.5 * ⟨z, H z⟩`.

Given the FW update `x_new = x - γ d`, this returns the exact minimizer
`γ* = ⟨∇f(x), d⟩ / ⟨d, H d⟩`, clamped to `[0, gamma_max]`.

Implementation detail:
- For homogeneous quadratics, `⟨d, H d⟩ = 2 f(d)`, so we can compute the curvature with one call to `f(d)`.

If you use this with a non-quadratic objective, the returned `γ` is not guaranteed to be optimal.
"""
struct QuadraticExactLineSearch <: FrankWolfe.LineSearchMethod end

FrankWolfe.build_linesearch_workspace(::QuadraticExactLineSearch, x, gradient) = nothing

function FrankWolfe.perform_line_search(
    ::QuadraticExactLineSearch,
    _,
    f,
    grad!,
    gradient,
    x,
    d,
    gamma_max,
    workspace,
    memory_mode,
)
    num = dot(gradient, d)
    if num <= 0
        return zero(num)
    end
    denom = 2 * f(d)
    if denom <= eps(float(denom))
        return zero(num)
    end
    gamma = num / denom
    if !isfinite(gamma)
        return zero(num)
    end
    return clamp(gamma, zero(gamma), gamma_max)
end

Base.print(io::IO, ::QuadraticExactLineSearch) = print(io, "QuadraticExact")