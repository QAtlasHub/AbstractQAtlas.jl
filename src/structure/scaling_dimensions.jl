# structure/scaling_dimensions.jl: the RG origin of the critical exponents,
# for the two kinds of fixed point this package can parameterise.
#
# Both follow one rule: DECLARE ONCE the independent data plus the spatial
# dimension `d`, DERIVE EVERYTHING else.  `ScalingDimensions` is the clean
# (or conventional random) case, two relevant eigenvalues; `InfiniteRandomness`
# is the activated case, where no finite `z` exists and the defining exponent
# is `ψ` instead.  `d` is a FIELD of each, not a per-call argument, because it
# is the one number a caller can supply wrongly without the answer looking
# wrong (see the warning on `InfiniteRandomness`).
#
# The four scaling laws in `relations/scaling.jl` (Rushbrooke, Widom,
# Fisher, Josephson) are written there as CHECKABLE identities among a
# supplied exponent set.  But they are not four independent physical
# facts: they are ALL consequences of a single structural statement — that
# the singular part of the free-energy density is a generalized
# homogeneous function of the reduced temperature `t` and the ordering
# field `h`,
#
#     f_s(t, h) = b^{-d} · f_s(b^{y_t} t, b^{y_h} h)                (∗)
#
# with just TWO relevant RG eigenvalues `y_t`, `y_h` and the spatial
# dimension `d`.  Every equilibrium exponent is a fixed rational function
# of `(y_t, y_h, d)`:
#
#     ν = 1/y_t,  α = 2 − d/y_t,  β = (d − y_h)/y_t,
#     γ = (2y_h − d)/y_t,  δ = y_h/(d − y_h),  η = d + 2 − 2y_h.
#
# Substituting these into the four laws, each collapses to `0` IDENTICALLY
# in `(y_t, y_h, d)` — so `critical_exponents(ScalingDimensions(...))`
# produces, by construction, a set that passes every scaling relation
# exactly.  This is the structural root the axioms are derived from: DECLARE
# ONCE (two eigenvalues + `d`), DERIVE EVERYTHING (the six exponents, hence
# the four laws).  No exponent VALUE is stored here — the eigenvalues are an
# input (measured, or read off an RG fixed point); the package owns only the
# universal `(y_t, y_h, d) ↦ exponents` map.
#
# Hyperscaling caveat: equation (∗) carries the bare `b^{-d}`, so the derived
# set satisfies the hyperscaling law (Josephson `2 − α = dν`, and the `d`
# in `η`) by construction.  Above the upper critical dimension a dangerous
# irrelevant variable spoils (∗) and the true (mean-field) exponents cease
# to obey Josephson except exactly at `d = d_upper`; there the eigenvalue
# parameterization no longer applies (see [`Josephson`](@ref)).

"""
    ScalingDimensions(y_t, y_h, d)

The renormalization-group data of a continuous transition: the thermal and
magnetic relevant eigenvalues `y_t`, `y_h` (the RG-flow exponents of the
reduced temperature and the ordering field) and the spatial dimension `d` of
the system this describes ([`SpatialDimension`](@ref)).  For a quantum critical
point that system is the classical image, so a chain is described here by its
2D image and `d = 2`; the chain's own `d = 1` is what
[`InfiniteRandomness`](@ref) and [`HarrisCriterion`](@ref) take.
These are the two-and-a-bit numbers the whole equilibrium exponent set is a
function of, via the homogeneity of the singular free energy
`f_s(t,h) = b^{-d} f_s(b^{y_t}t, b^{y_h}h)` — see
[`critical_exponents`](@ref).

Arguments are promoted to a common type; pass `Rational`s for an exact set
(`ScalingDimensions(1//1, 15//8, 2)` is 2D Ising) or `Float64` for a
numerical fixed point.

```julia
critical_exponents(ScalingDimensions(1//1, 15//8, 2))
# (α = 0//1, β = 1//8, γ = 7//4, δ = 15//1, ν = 1//1, η = 1//4)   ← 2D Ising, exact
```
"""
struct ScalingDimensions{T<:Real}
    y_t::T
    y_h::T
    d::T
end
function ScalingDimensions(y_t::Real, y_h::Real, d::Real)
    return ScalingDimensions(promote(y_t, y_h, d)...)
end
export ScalingDimensions

"""
    critical_exponents(s::ScalingDimensions) -> NamedTuple

The full equilibrium exponent set `(α, β, γ, δ, ν, η)` DERIVED from the
RG eigenvalues in `s` — nothing hand-entered:

- `ν = 1/y_t`                (correlation length, `ξ ∼ |t|^{-ν}`)
- `α = 2 − d/y_t`            (specific heat, `C ∼ |t|^{-α}`)
- `β = (d − y_h)/y_t`        (order parameter, `M ∼ |t|^{+β}`)
- `γ = (2y_h − d)/y_t`       (susceptibility, `χ ∼ |t|^{-γ}`)
- `δ = y_h/(d − y_h)`        (critical isotherm, `M ∼ h^{1/δ}`)
- `η = d + 2 − 2y_h`         (anomalous dimension, `G(r) ∼ r^{-(d-2+η)}`)

With `Rational` eigenvalues the result is exact, and it satisfies
[`Rushbrooke`](@ref), [`Widom`](@ref), [`Fisher`](@ref) and
[`Josephson`](@ref) with residual `≡ 0` for *any* `s`
([`exponents_consistent`](@ref) is `true` by construction).
"""
function critical_exponents(s::ScalingDimensions)
    yt, yh, d = s.y_t, s.y_h, s.d
    return (
        α=2 - d / yt,
        β=(d - yh) / yt,
        γ=(2yh - d) / yt,
        δ=yh / (d - yh),
        ν=1 / yt,
        η=d + 2 - 2yh,
    )
end
export critical_exponents

"""
    critical_exponent(name::Symbol, s::ScalingDimensions) -> Real

A single derived exponent (`:α`, `:β`, `:γ`, `:δ`, `:ν`, or `:η`) of `s`
— a keyed view of [`critical_exponents`](@ref).
"""
critical_exponent(name::Symbol, s::ScalingDimensions) = critical_exponents(s)[name]
export critical_exponent

"""
    scaling_dimensions(; ν, η, d) -> ScalingDimensions

Invert the exponent map: recover the RG eigenvalues from the two
independent exponents that fix them, `y_t = 1/ν` and `y_h = (d + 2 − η)/2`,
at dimension `d`.  Composing with [`critical_exponents`](@ref) closes the
loop — every other exponent (`α, β, γ, δ`) is then reconstructed from just
`(ν, η, d)`, a direct expression of the two-eigenvalue structure:

```julia
s = scaling_dimensions(ν = 1//1, η = 1//4, d = 2)   # 2D Ising eigenvalues
critical_exponents(s).δ                              # 15//1  (δ from ν, η, d alone)
```
"""
function scaling_dimensions(; ν, η, d)
    y_t = 1 / ν
    y_h = (d + 2 - η) / 2
    return ScalingDimensions(y_t, y_h, d)
end
export scaling_dimensions

# ─── Infinite-randomness fixed points ────────────────────────────────────
#
# At an infinite-randomness fixed point the free energy is not a homogeneous
# function of `(t, h)` with two eigenvalues: the dynamics is activated, no
# finite `z` exists, and the scale that plays `y_t`'s role is `ψ`.  What
# survives is the same DECLARE-ONCE structure with a different independent
# set, `(ψ, ν, x_m)` plus `d`, from which the rest of the exponent table
# follows through relations this package already states:
#
#     β     = ν·x_m                (OrderParameterDimension,   Eq. (A.11))
#     ν_typ = ν·(1 − ψ)            (TypicalCorrelationLength,  Eq. (9.4))
#     φ     = (d − x_m)/ψ          (ActivatedMomentGrowth,     Eq. (A.21))
#
# Surface exponents are NOT derivable from the bulk set, here as in the clean
# case: `x_m^s` is independent data and is not carried.

"""
    InfiniteRandomness(ψ, ν, x_m, d)

The data of an infinite-randomness fixed point: the activated exponent `ψ`
(`ln(1/Δ) ∼ ξ^ψ`), the average correlation-length exponent `ν`, the bulk
order-parameter scaling dimension `x_m`, and the SPATIAL dimension `d`.  The
counterpart of [`ScalingDimensions`](@ref) for a fixed point that has no finite
dynamical exponent, and the same contract: these are the inputs, every other
exponent is derived by [`critical_exponents`](@ref).

!!! warning "`d` is the chain's own dimension, not its classical image's"
    A random transverse-field Ising chain is `d = 1` here.  Its clean critical
    point maps to the 2D classical Ising model, and atlases hand out `d = 2`
    for that table, so 2 is the number nearest to hand and it is the wrong one:
    `φ` comes out 3.618 instead of the golden mean 1.618, and nothing flags it.
    Quenched disorder is constant along imaginary time, so it lives in the
    chain's `d`, not the image's `d + z` (see [`SpatialDimension`](@ref)).
    Holding it in the struct makes that one decision at construction rather
    than one per call.  There is no `d + z` to confuse it with anyway: `z` is
    infinite at an infinite-randomness fixed point.

Arguments are promoted to a common type; pass `Rational`s where the values are
rational.  `ψ ≤ 0` is refused rather than accepted: it names a conventional
fixed point, which is [`ScalingDimensions`](@ref)'s job, and it would divide by
zero in `φ`.

```julia
critical_exponents(InfiniteRandomness(1//2, 2//1, (3-sqrt(5))/4, 1))
# (β = 0.381…, ν = 2.0, ν_typ = 1.0, ψ = 0.5, x_m = 0.190…, φ = 1.618…)
```
"""
struct InfiniteRandomness{T<:Real}
    ψ::T
    ν::T
    x_m::T
    d::T
    function InfiniteRandomness(ψ::T, ν::T, x_m::T, d::T) where {T<:Real}
        ψ > 0 || error(
            "InfiniteRandomness: ψ = $ψ is not an infinite-randomness fixed point. " *
            "ψ > 0 is what the name means; a fixed point with ψ = 0 has a finite " *
            "dynamical exponent and is parameterised by ScalingDimensions.",
        )
        ν > 0 || error("InfiniteRandomness: ν = $ν must be positive.")
        d > 0 || error("InfiniteRandomness: d = $d must be positive.")
        return new{T}(ψ, ν, x_m, d)
    end
end
function InfiniteRandomness(ψ::Real, ν::Real, x_m::Real, d::Real)
    return InfiniteRandomness(promote(ψ, ν, x_m, d)...)
end
export InfiniteRandomness

"""
    critical_exponents(s::InfiniteRandomness) -> NamedTuple

The exponent set `(β, ν, ν_typ, ψ, x_m, φ)` of an infinite-randomness fixed
point, DERIVED from `s`: the three independent inputs plus `d`, nothing
hand-entered:

- `β     = ν·x_m`         (order parameter, `m ∼ |δ|^{+β}`)
- `ν_typ = ν·(1 − ψ)`     (typical correlation length, always `< ν` for `ψ > 0`)
- `φ     = (d − x_m)/ψ`   (cluster-moment growth, `μ ∼ |ln Ω|^φ`)

`ψ` and `x_m` are inputs that are themselves exponents, so they appear in the
set; `d` does not, matching [`critical_exponents`](@ref)`(::ScalingDimensions)`.
Use [`exponents_consistent`](@ref)`(s)` to sweep the registry without having to
restate it.

By construction the result satisfies [`OrderParameterDimension`](@ref),
[`TypicalCorrelationLength`](@ref) and [`ActivatedMomentGrowth`](@ref) with
residual exactly zero, for any `s`.
"""
function critical_exponents(s::InfiniteRandomness)
    ψ, ν, x_m, d = s.ψ, s.ν, s.x_m, s.d
    return (β=ν * x_m, ν=ν, ν_typ=ν * (1 - ψ), ψ=ψ, x_m=x_m, φ=(d - x_m) / ψ)
end

"""
    critical_exponent(name::Symbol, s::InfiniteRandomness) -> Real

A single derived exponent (`:β`, `:ν`, `:ν_typ`, `:ψ`, `:x_m` or `:φ`) of `s`.
"""
critical_exponent(name::Symbol, s::InfiniteRandomness) = critical_exponents(s)[name]

"""
    infinite_randomness(; ν, ν_typ, x_m, d) -> InfiniteRandomness

Invert the exponent map: recover the defining `ψ = 1 − ν_typ/ν` from the two
correlation lengths, at dimension `d`.  Composing with
[`critical_exponents`](@ref) closes the loop, the way
[`scaling_dimensions`](@ref) does for the clean case.

```julia
s = infinite_randomness(; ν = 2//1, ν_typ = 1//1, x_m = 1//4, d = 1)
s.ψ            # 1//2
```
"""
function infinite_randomness(; ν, ν_typ, x_m, d)
    return InfiniteRandomness(1 - ν_typ / ν, ν, x_m, d)
end
export infinite_randomness
