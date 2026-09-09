# ext/AbstractQAtlasForwardDiffExt.jl — AD realization of the response
# genealogy via ForwardDiff.
#
# Each method evaluates the derivative a supplied-derivative relation
# needs directly from the potential function, so a downstream calculation
# can go from "a free-energy function F(h)" to "the magnetization" (and
# on to `check(MagnetizationResponse(); …)`) without hand-coding the
# derivative.

module AbstractQAtlasForwardDiffExt

using AbstractQAtlas
import AbstractQAtlas: peierls_current, thermal_derivative   # extended below → must import
using AbstractQAtlas:
    response_order,
    indices,
    derivative_edge,
    derivative_order,
    potential_root,
    FreeEnergy,
    GrandPotential,
    Susceptibility,
    SpecificHeat,
    Energy
using ForwardDiff: derivative

# n-th derivative of a scalar function by nested ForwardDiff (n small —
# response orders are 1–3).
_nth(f, x, n::Integer) = n == 0 ? f(x) : _nth(y -> derivative(f, y), x, n - 1)

# ── generic, genealogy-driven response ────────────────────────────────────
# For a FreeEnergy-rooted response that is a pure single-field derivative of
# the free energy — the extensive-variable tree `M = −∂F/∂h`, `S = −∂F/∂T`, and
# any higher single-field response — read the derivative ORDER and FIELD from
# the genealogy (`derivative_edge` / `derivative_order`) instead of hand-coding
# each quantity.  A new such quantity is AD-covered the moment it declares its
# edge — the map is not hand-maintained (core-functions.md, pillar 4).  The net
# sign is −1 across this tree (the `F → response` edge is `−∂/∂field`, the
# intra-response edges are `+∂/∂field`).  Irregular members of the genealogy
# (`Energy`, whose potential is `βF`, and `SpecificHeat`, reached through `U`)
# keep the explicit methods below; `Susceptibility` keeps its own method for the
# off-diagonal guard.  `x::Number` on every point argument keeps this generic
# method and the concrete ones unambiguous.
function thermal_derivative(q::AbstractQuantity, F, x::Number)
    e = derivative_edge(q)
    e === nothing && error(
        "thermal_derivative: $(typeof(q)) is not a response function " *
        "(no derivative_edge) — nothing to differentiate.",
    )
    root = potential_root(q)
    (root === FreeEnergy || root === GrandPotential) || error(
        "thermal_derivative: the generic genealogy path handles the standard response " *
        "potentials (FreeEnergy `M=−∂F/∂h`, GrandPotential `N=−∂Ω/∂μ`); $(typeof(q)) " *
        "roots at $root — use its explicit method (e.g. Energy/SpecificHeat).",
    )
    return -_nth(F, x, derivative_order(q, e.field()))   # e.field is a TYPE → instantiate
end

# χ⁽ⁿ⁾_{α;β₁…βₙ} = −∂ⁿ⁺¹F/∂h_α∂h_{β₁}…∂h_{βₙ}.  With a SINGLE-field function
# F(h) only the DIAGONAL component (all indices equal) is defined — an
# off-diagonal component needs partials w.r.t. distinct field directions,
# so guard against silently returning the diagonal for an off-diagonal ask.
function thermal_derivative(χ::Susceptibility, F, h::Number)
    idx = indices(χ)
    all(==(idx[1]), idx) || error(
        "thermal_derivative(::Susceptibility, F, h::Number) with a single-field " *
        "function computes only the DIAGONAL χ⁽ⁿ⁾ (all indices equal); got " *
        "off-diagonal $(idx). Pass a multi-field potential F(h⃗) and the field-" *
        "component ordering: thermal_derivative(χ, F, h⃗, components).",
    )
    return -_nth(F, h, response_order(χ) + 1)
end

# Multi-field / off-diagonal: F is a function of a field VECTOR `h⃗`, and
# `components[k]` names the field direction of slot k.  The tensor
# susceptibility is the mixed partial over the directions in `indices(χ)`
# (response index α included: M_α = −∂F/∂h_α):
#   χ⁽ⁿ⁾_{α;β₁…βₙ} = −∂ⁿ⁺¹F / ∂h_α ∂h_{β₁} … ∂h_{βₙ}
# with the diagonal (all-equal indices) reproducing the single-field result.
function _slot(s, components)
    p = findfirst(==(s), components)
    p === nothing && error(
        "thermal_derivative: field component $(repr(s)) not in components $(components)"
    )
    return p
end

_bump(h⃗, p, t) = [i == p ? h⃗[i] + t : h⃗[i] for i in eachindex(h⃗)]

_partial(F, h⃗, ::Tuple{}) = F(h⃗)
function _partial(F, h⃗, slots::Tuple)
    p = first(slots)
    return derivative(t -> _partial(F, _bump(h⃗, p, t), Base.tail(slots)), zero(eltype(h⃗)))
end

function thermal_derivative(χ::Susceptibility, F, h⃗::AbstractVector, components)
    slots = map(s -> _slot(s, components), indices(χ))
    return -_partial(F, h⃗, slots)
end

# S = −∂F/∂T is handled by the generic genealogy method above (order-1 single-
# field response of the free energy w.r.t. Temperature).

# ── irregular members: not a pure single-field power of F ──────────────────
# C = ∂U/∂T  (reached through U, not F directly)
thermal_derivative(::SpecificHeat, U, T::Number) = derivative(U, T)

# U = ∂(βF)/∂β  (Gibbs–Helmholtz; the potential is βF, not F — pass βF(β))
thermal_derivative(::Energy, βF, β::Number) = derivative(βF, β)

# J = -∂H/∂A. The same forward-mode route `thermal_derivative` takes, and the same sign
# convention as `M = -∂F/∂h`: differentiating at the seam is what keeps a consumer from
# hand-coding a current operator and getting its sign wrong in a way no static check sees.
peierls_current(H_of_A, A::Number) = -derivative(H_of_A, A)

end # module AbstractQAtlasForwardDiffExt
