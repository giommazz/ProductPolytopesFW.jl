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

