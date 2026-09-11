# structure/relevance.jl — does a perturbation change the fixed point?
#
# A criterion is neither a relation nor a bound.  A relation's residual must
# vanish; a bound's must keep a sign.  Here BOTH signs are answers, and the
# zero is a third one — so these deliberately do not enter the relation
# registry, where `check_all` would report "failed" for every system that is
# simply not marginal.
#
# References: Harris, [Harris1974](@cite); Luck, [Luck1993](@cite);
# Weinrib–Halperin, [WeinribHalperin1983](@cite).  The three are collected,
# with the equation numbers cited below, in Iglói–Monthus,
# [IgloiMonthus2005](@cite) (arXiv `cond-mat/0502448`).

"""
    RelevanceCriterion

A predicate on critical exponents deciding whether a perturbation changes the
fixed point.  Answered by [`relevance`](@ref); the underlying signed quantity is
[`margin`](@ref).
"""
abstract type RelevanceCriterion end
export RelevanceCriterion

"""
    margin(c::RelevanceCriterion; kwargs...) -> Real

The signed distance from marginality, in the convention `> 0` irrelevant,
`0` marginal, `< 0` relevant.  Same sign convention for every criterion, so
`relevance` needs no per-criterion orientation.
"""
function margin end
export margin

"""
    relevance(c::RelevanceCriterion; atol=0, kwargs...) -> Symbol

`:irrelevant`, `:marginal` or `:relevant`, from the sign of [`margin`](@ref).

`atol` widens the marginal band; it defaults to `0`, so an exponent set that is
marginal in exact arithmetic reports `:marginal` and a floating-point one very
likely will not.  Pass a tolerance when the exponents are estimates.
"""
function relevance(c::RelevanceCriterion; atol::Real=0, kwargs...)
    m = margin(c; kwargs...)
    abs(m) <= atol && return :marginal
    return m > 0 ? :irrelevant : :relevant
end
export relevance

"""
    HarrisCriterion() <: RelevanceCriterion

Uncorrelated quenched disorder is irrelevant at a clean fixed point with
correlation-length exponent `ν₀` in `d` dimensions when

`ν₀ > 2/d`,

so `margin = ν₀ − 2/d`.  `ν₀` is the exponent of the PURE system: the criterion
asks what the clean fixed point does to a perturbation, not what the disordered
one looks like.

The clean transverse-field Ising chain has `ν₀ = 1` at `d = 1`, giving
`margin = −1`: disorder is relevant, which is why the random chain flows to an
infinite-randomness fixed point instead.

Reference: [Harris1974](@cite); stated as Eq. (5.9) of [IgloiMonthus2005](@cite).
"""
struct HarrisCriterion <: RelevanceCriterion end
export HarrisCriterion

margin(::HarrisCriterion; ν₀, d, _extra...) = ν₀ - 2 / d

"""
    LuckCriterion() <: RelevanceCriterion

Luck's extension to APERIODIC modulation, where the perturbation is deterministic
and its strength is set by how its fluctuations grow, `Δ(L) ∼ L^ω`.  Irrelevant
when

`ν₀ > 1/(1 − ω)`,

so `margin = ν₀ − 1/(1 − ω)`.  `ω = 1/2` is a random sequence and `ω = −1` the
Fibonacci one, whose bounded fluctuations make it irrelevant wherever `ν₀ > 1/2`.

At `ω = 1/2` this IS [`HarrisCriterion`](@ref) in one dimension — `1/(1−1/2) = 2
= 2/d` — which is the consistency the two must have and a way to check either.

Reference: [Luck1993](@cite); stated as Eq. (10.9) of [IgloiMonthus2005](@cite),
with the wandering exponent defined in its Eq. (10.8).
"""
struct LuckCriterion <: RelevanceCriterion end
export LuckCriterion

margin(::LuckCriterion; ν₀, ω, _extra...) = ν₀ - 1 / (1 - ω)

"""
    WeinribHalperinCriterion() <: RelevanceCriterion

Disorder whose correlator decays as a power, `G(r) ∼ r^{−ρ}`, rather than being
uncorrelated.  Correlations are irrelevant — the uncorrelated universality class
survives — when

`ρ > 2/ν`,

so `margin = ρ − 2/ν`.  Here `ν` is the correlation-length exponent of the
UNCORRELATED disordered fixed point, not the clean one that
[`HarrisCriterion`](@ref) takes: this criterion asks a later question, whether
correlations move a fixed point disorder has already changed.

Reference: [WeinribHalperin1983](@cite); stated as Eq. (10.2) of
[IgloiMonthus2005](@cite).
"""
struct WeinribHalperinCriterion <: RelevanceCriterion end
export WeinribHalperinCriterion

margin(::WeinribHalperinCriterion; ν, ρ, _extra...) = ρ - 2 / ν
