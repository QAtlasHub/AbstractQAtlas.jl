# structure/relevance.jl — does a perturbation change the fixed point?
#
# A criterion is neither a relation nor a bound.  A relation's residual must
# vanish; a bound's must keep a sign.  Here BOTH signs are answers, and the
# zero is a third one — so these deliberately do not enter the relation
# registry, where `check_all` would report "failed" for every system that is
# simply not marginal.  Measured: as a `@bound`, `check(...; ν₀=1, d=1)` — the
# clean Ising chain, a textbook RELEVANT case — comes back `false`, and
# `check_all` with it.
#
# The cost, which is real: staying outside `AbstractRelation` also forfeits
# `domain`, the `Bag`/`VariableKey` front door and load-time validation.  The
# type-keyed bag is this package's structural answer to "two exponents that look
# alike", and these criteria carry exactly such a pair — hence `ν_dis` below.
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

# Without this a criterion missing its `margin` surfaces as a bare MethodError at
# first use — possibly downstream, far from where it was declared. Same shape as
# `critical_scaling`'s default-then-named-check next door in `criticality.jl`.
function margin(c::RelevanceCriterion; kwargs...)
    return error(
        "margin: $(nameof(typeof(c))) is a `RelevanceCriterion` with no `margin` " *
        "method. Define one returning the signed distance from marginality — " *
        "positive irrelevant, zero marginal, negative relevant.",
    )
end

"""
    relevance(c::RelevanceCriterion; atol=0, kwargs...) -> Symbol

`:irrelevant`, `:marginal` or `:relevant`, from the sign of [`margin`](@ref).

`atol` widens the marginal band; it defaults to `0`, so an exponent set that is
marginal in exact arithmetic reports `:marginal` and a floating-point one very
likely will not.  Pass a tolerance when the exponents are estimates.

`atol` is reserved: it is consumed here and never reaches [`margin`](@ref), so a
criterion may not name a physics variable `atol`.
"""
function relevance(c::RelevanceCriterion; atol::Real=0, kwargs...)
    m = margin(c; kwargs...)
    # NaN would otherwise be reported as `:relevant`: IEEE makes BOTH `abs(m) <= atol`
    # and `m > 0` false, so the ternary's else-branch wins and an upstream failure
    # becomes a confident physical claim. No `atol` rescues it, not even `Inf`.
    isfinite(m) ||
        error("relevance: $(nameof(typeof(c))) gave a non-finite margin ($m). That is \
               an unusable input reaching a criterion, not a verdict — check the \
               exponents rather than reading a relevance from it.")
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

function margin(::HarrisCriterion; ν₀, d, _extra...)
    # `d > 0` is not pedantry: at d = 0 the margin is ±Inf, and IEEE signed zero
    # makes `d = -0.0` and `d = +0.0` give OPPOSITE verdicts for the same
    # "zero dimensions".
    d > 0 || throw(ArgumentError("HarrisCriterion: d must be > 0; got $d"))
    return ν₀ - 2 / d
end

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

function margin(::LuckCriterion; ν₀, ω, _extra...)
    # ω = 1 divides by zero, and ω > 1 means fluctuations outgrowing the system,
    # which the criterion is not derived for — random is 1/2 and Fibonacci is −1.
    ω < 1 || throw(
        ArgumentError(
            "LuckCriterion: the wandering exponent must be < 1; got $ω. At ω = 1 the " *
            "bound diverges and above it the fluctuations outgrow L, which this " *
            "criterion does not describe.",
        ),
    )
    return ν₀ - 1 / (1 - ω)
end

"""
    WeinribHalperinCriterion() <: RelevanceCriterion

Disorder whose correlator decays as a power, `G(r) ∼ r^{−ρ}`, rather than being
uncorrelated.  Correlations are irrelevant — the uncorrelated universality class
survives — when

`ρ > 2/ν`,

so `margin = ρ − 2/ν_dis`.  The keyword is `ν_dis`, not `ν`, on purpose: this is
the correlation-length exponent of the UNCORRELATED DISORDERED fixed point, one
subscript away from [`HarrisCriterion`](@ref)'s clean `ν₀` and a different
number.  Passing a clean exponent here is syntactically fine and physically
wrong, so the two are not allowed to look alike.

The criteria ask different questions in sequence: Harris, whether disorder
matters at all; this one, whether correlations move a fixed point disorder has
already changed.

Reference: [WeinribHalperin1983](@cite); stated as Eq. (10.2) of
[IgloiMonthus2005](@cite).
"""
struct WeinribHalperinCriterion <: RelevanceCriterion end
export WeinribHalperinCriterion

function margin(::WeinribHalperinCriterion; ν_dis, ρ, _extra...)
    ν_dis > 0 ||
        throw(ArgumentError("WeinribHalperinCriterion: ν_dis must be > 0; got $ν_dis"))
    ρ > 0 || throw(
        ArgumentError("WeinribHalperinCriterion: the decay exponent ρ must be > 0; got $ρ"),
    )
    return ρ - 2 / ν_dis
end
