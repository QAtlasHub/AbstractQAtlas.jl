# structure/transitions.jl — phase-transition classification.
#
# The generic, model-independent facts about *kinds* of transition:
# what distinguishes first-order from continuous from BKT.  The
# Ehrenfest picture roots the distinction in the free energy — a
# transition is n-th order when the n-th derivative of the free energy
# is the first one to be discontinuous / singular.  Every other trait
# (order parameter? latent heat? critical exponents?) follows from that
# one structural fact and is declared here so downstream code can
# dispatch on transition kind instead of hard-coding case analysis.

"""
    AbstractTransition

Abstract parent type for phase-transition kinds.  Concrete subtypes
([`FirstOrder`](@ref), [`ContinuousTransition`](@ref),
[`KosterlitzThouless`](@ref)) are singleton tags carrying the generic
classification traits ([`ehrenfest_order`](@ref),
[`has_order_parameter`](@ref), [`has_latent_heat`](@ref),
[`has_critical_exponents`](@ref)).

The classification is anchored in the free energy: an `n`-th order
transition is one whose `n`-th free-energy derivative is the first to be
discontinuous or divergent.
"""
abstract type AbstractTransition end
export AbstractTransition

"""
    FirstOrder <: AbstractTransition

A first-order (discontinuous) transition: the first derivative of the
free energy jumps.  Concretely — a discontinuity in the order parameter
`M = −∂F/∂h` and/or the entropy `S = −∂F/∂T` (the latter giving a latent
heat `L = T ΔS`).  No diverging correlation length, hence no critical
exponents in the continuous-transition sense.
"""
struct FirstOrder <: AbstractTransition end
export FirstOrder

"""
    ContinuousTransition <: AbstractTransition

A continuous (second-order / critical) transition: the free energy and
its first derivatives are continuous, while second derivatives — the
[`SpecificHeat`](@ref) `C = −T ∂²F/∂T²` and the
susceptibility `χ = −∂²F/∂h²` — diverge.  The correlation length
diverges, so the singularities are power laws governed by
[`CriticalExponents`](@ref); the [`critical_scaling`](@ref)
correspondence assigns each observable its exponent.
"""
struct ContinuousTransition <: AbstractTransition end
export ContinuousTransition

"""
    KosterlitzThouless <: AbstractTransition

The Berezinskii–Kosterlitz–Thouless transition (2D XY and relatives):
an infinite-order transition with an essential singularity in the free
energy (`ξ ∼ exp(c/√(T−T_c))`, not a power law) and **no local order
parameter** (Mermin–Wagner).  The standard equilibrium critical
exponents do not apply; the transition is characterized instead by the
universal helicity-modulus jump.
"""
struct KosterlitzThouless <: AbstractTransition end
export KosterlitzThouless

# ─── Classification traits ──────────────────────────────────────────────

"""
    ehrenfest_order(::AbstractTransition) -> Union{Int,Float64}

The order of the transition in the (generalized) Ehrenfest sense: the
index of the lowest free-energy derivative that is discontinuous or
singular.  `1` for [`FirstOrder`](@ref), `2` for
[`ContinuousTransition`](@ref), `Inf` for the essential singularity of
[`KosterlitzThouless`](@ref).
"""
ehrenfest_order(::FirstOrder) = 1
ehrenfest_order(::ContinuousTransition) = 2
ehrenfest_order(::KosterlitzThouless) = Inf
export ehrenfest_order

"""
    has_order_parameter(::AbstractTransition) -> Bool

Whether the transition is characterized by a local order parameter that
is nonzero on one side.  `false` for [`KosterlitzThouless`](@ref) (a
topological transition with only quasi-long-range order).
"""
has_order_parameter(::FirstOrder) = true
has_order_parameter(::ContinuousTransition) = true
has_order_parameter(::KosterlitzThouless) = false
export has_order_parameter

"""
    has_latent_heat(::AbstractTransition) -> Bool

Whether the transition releases a latent heat `L = T ΔS` — i.e. whether
the entropy `S = −∂F/∂T` is discontinuous.  Only [`FirstOrder`](@ref).
"""
has_latent_heat(::FirstOrder) = true
has_latent_heat(::ContinuousTransition) = false
has_latent_heat(::KosterlitzThouless) = false
export has_latent_heat

"""
    has_critical_exponents(::AbstractTransition) -> Bool

Whether power-law critical exponents in the [`CriticalExponents`](@ref)
sense apply — true only for [`ContinuousTransition`](@ref), where the
correlation length diverges as a power law.
"""
has_critical_exponents(::FirstOrder) = false
has_critical_exponents(::ContinuousTransition) = true
has_critical_exponents(::KosterlitzThouless) = false
export has_critical_exponents
