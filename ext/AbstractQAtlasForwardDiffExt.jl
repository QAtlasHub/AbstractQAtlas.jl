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
import AbstractQAtlas: thermal_derivative        # extended below → must import
using AbstractQAtlas:
    response_order,
    indices,
    Magnetization,
    Susceptibility,
    ThermalEntropy,
    SpecificHeat,
    Energy
using ForwardDiff: derivative

# n-th derivative of a scalar function by nested ForwardDiff (n small —
# response orders are 1–3).
_nth(f, x, n::Integer) = n == 0 ? f(x) : _nth(y -> derivative(f, y), x, n - 1)

# M_α = −∂F/∂h  (first field-derivative of the free energy)
thermal_derivative(::Magnetization, F, h) = -derivative(F, h)

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

# S = −∂F/∂T
thermal_derivative(::ThermalEntropy, F, T) = -derivative(F, T)

# C = ∂U/∂T
thermal_derivative(::SpecificHeat, U, T) = derivative(U, T)

# U = ∂(βF)/∂β  (Gibbs–Helmholtz; pass the βF function of β)
thermal_derivative(::Energy, βF, β) = derivative(βF, β)

end # module AbstractQAtlasForwardDiffExt
