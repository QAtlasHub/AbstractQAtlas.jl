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
import AbstractQAtlas: nth_derivative, peierls_current, thermal_derivative
using AbstractQAtlas:
    _genealogy_derivative,
    _susceptibility_derivative,
    indices,
    Susceptibility,
    SpecificHeat,
    Energy
using ForwardDiff: derivative

# n-th derivative of a scalar function by nested ForwardDiff (n small —
# response orders are 1–3).
_nth(f, x, n::Integer) = n == 0 ? f(x) : _nth(y -> derivative(f, y), x, n - 1)

# `AutoDiff` is a route like the others: this is the method whose absence made
# every route-taking entry point special-case it by name.
nth_derivative(::AbstractQAtlas.AutoDiff, f, x, n::Integer) = _nth(f, x, n)

# ── generic, genealogy-driven response ────────────────────────────────────
# Which order, which field and the net sign all come from `_genealogy_derivative`
# (src/derivative_routes.jl), shared with the finite-difference routes so the two
# cannot drift. This supplies the AD `nth` and nothing else. `x::Number` keeps
# this generic method and the concrete ones below unambiguous.
function thermal_derivative(q::AbstractQuantity, F, x::Number)
    return _genealogy_derivative(q, F, x, _nth)
end

# χ⁽ⁿ⁾_{α;β₁…βₙ} = −∂ⁿ⁺¹F/∂h_α∂h_{β₁}…∂h_{βₙ}, diagonal only from a single-field
# F(h). The guard and the order live in `_susceptibility_derivative`, shared with
# the finite-difference routes; this supplies the AD `nth`.
function thermal_derivative(χ::Susceptibility, F, h::Number)
    return _susceptibility_derivative(χ, F, h, _nth)
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
function peierls_current(H_of_A, A::AbstractQAtlas.VectorPotential{1})
    return -derivative(a -> H_of_A(AbstractQAtlas.VectorPotential(a)), A[1])
end

end # module AbstractQAtlasForwardDiffExt
