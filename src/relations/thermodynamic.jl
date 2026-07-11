# relations/thermodynamic.jl — fluctuation–dissipation identities.
#
# Conventions match the QAtlas identities plane (test/identities/):
#
#   SpecificHeat     c_v = β² · Var(E) / N        (E = TOTAL energy)
#   Susceptibility   χ   = β  · Var(M) / N        (M = TOTAL magnetization)
#
# i.e. both responses are PER SITE when `N` is the number of sites and
# the fluctuating quantity is extensive.  Pass `N = 1` (the default) to
# work with total responses, or with already-intensive variances.
#
# `solve(rel, Val(:C); ...)` / `solve(rel, Val(:χ); ...)` double as the
# estimators downstream Monte-Carlo / ED packages should use, so the
# formula lives in exactly one place.

"""
    SpecificHeatFDT <: AbstractRelation

The energy fluctuation–dissipation identity

`c_v = β² (⟨E²⟩ − ⟨E⟩²) / N = β² Var(E) / N`,

with `E` the total energy and `N` the site count (`N = 1` ⇒ total
specific heat).  Equivalently `c_v = −β² ∂⟨E⟩/∂β / N`: the fluctuation
route and the temperature-response route must agree — that is the
content of the relation, and how it is tested.

Variables: `C`, `var_E`, `β` (or `T`), `N` (default 1).

```julia
residual(SpecificHeatFDT(); C=c, var_E=v, β=β, N=N)   # c − β²v/N
solve(SpecificHeatFDT(), Val(:C); var_E=v, β=β, N=N)  # the estimator
```
"""
struct SpecificHeatFDT <: AbstractRelation end
export SpecificHeatFDT

_beta(; β=nothing, T=nothing) =
    if β !== nothing
        β
    elseif T !== nothing
        1 / T
    else
        error("pass either β or T")
    end

function residual(::SpecificHeatFDT; C, var_E, β=nothing, T=nothing, N=1)
    b = _beta(; β=β, T=T)
    return C - b^2 * var_E / N
end
function solve(::SpecificHeatFDT, ::Val{:C}; var_E, β=nothing, T=nothing, N=1)
    b = _beta(; β=β, T=T)
    return b^2 * var_E / N
end
function solve(::SpecificHeatFDT, ::Val{:var_E}; C, β=nothing, T=nothing, N=1)
    b = _beta(; β=β, T=T)
    return C * N / b^2
end

"""
    SusceptibilityFDT <: AbstractRelation

The static (zero-frequency) fluctuation–dissipation identity

`χ = β (⟨M²⟩ − ⟨M⟩²) / N = β Var(M) / N`,

with `M` the total (extensive) order parameter conjugate to the applied
field and `N` the site count (`N = 1` ⇒ total susceptibility).  This is
the h → 0 limit of the linear response `χ = ∂⟨M⟩/∂h` for
`H(h) = H₀ − h·M`; response route and fluctuation route must agree.

Variables: `χ`, `var_M`, `β` (or `T`), `N` (default 1).
"""
struct SusceptibilityFDT <: AbstractRelation end
export SusceptibilityFDT

function residual(::SusceptibilityFDT; χ, var_M, β=nothing, T=nothing, N=1)
    b = _beta(; β=β, T=T)
    return χ - b * var_M / N
end
function solve(::SusceptibilityFDT, ::Val{:χ}; var_M, β=nothing, T=nothing, N=1)
    b = _beta(; β=β, T=T)
    return b * var_M / N
end
function solve(::SusceptibilityFDT, ::Val{:var_M}; χ, β=nothing, T=nothing, N=1)
    b = _beta(; β=β, T=T)
    return χ * N / b
end

"""
    LinearResponseFDT <: AbstractRelation

The general static linear-response identity for a perturbation
`H(λ) = H₀ − λ·O` with `[O, H₀] = 0` (classical statistics or a
commuting observable):

`∂⟨O⟩/∂λ = β Var(O)`.

Variables: `dO_dλ`, `var_O`, `β` (or `T`).
[`SusceptibilityFDT`](@ref) is this relation with `O = M` and a `1/N`
normalization.
"""
struct LinearResponseFDT <: AbstractRelation end
export LinearResponseFDT

function residual(::LinearResponseFDT; dO_dλ, var_O, β=nothing, T=nothing)
    b = _beta(; β=β, T=T)
    return dO_dλ - b * var_O
end
function solve(::LinearResponseFDT, ::Val{:dO_dλ}; var_O, β=nothing, T=nothing)
    b = _beta(; β=β, T=T)
    return b * var_O
end
