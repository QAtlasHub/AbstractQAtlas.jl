# relations/scaling.jl — critical-exponent scaling laws.
#
# Each law is one @relation declaration; residual/check/solve for every
# variable follow mechanically (all four laws are affine in every
# variable).  Exact exponent sets (2D Ising rationals, mean-field at the
# upper critical dimension) satisfy them with residual ≡ 0 in exact
# arithmetic; numerical sets (3D Ising bootstrap) within quoted errors.
#
# Each relation carries its own source.  Equation, table and section numbers
# for [IgloiMonthus2005](@cite) are those of the arXiv version
# (`cond-mat/0502448`) and were read off it, not counted from the LaTeX — that
# review numbers by section, so a count gives different numbers entirely.
# The four classical identities are cited at paper level: each is the whole
# point of a short paper, and none of them is open access to locate within.

"""
    Rushbrooke <: AbstractRelation

The Rushbrooke identity `α + 2β + γ = 2`.

Reference: [Rushbrooke1963](@cite) (J. Chem. Phys. **39**, 842), where it is
derived thermodynamically as the INEQUALITY `α′ + 2β + γ′ ≥ 2`; equality is the
scaling hypothesis, not the 1963 result.

```julia
residual(Rushbrooke(); α=0//1, β=1//8, γ=7//4)   # == 0//1 (2D Ising, exact)
solve(Rushbrooke(), Val(:γ); α=0//1, β=1//8)     # == 7//4
```
"""
@relation :scaling Rushbrooke(α, β, γ) = α + 2β + γ - 2

"""
    Widom <: AbstractRelation

The Widom identity `γ = β(δ − 1)`.

Reference: [Widom1965](@cite) (J. Chem. Phys. **43**, 3898), from the
homogeneous equation of state.
"""
@relation :scaling Widom(β, γ, δ) = γ - β * (δ - 1)

"""
    Fisher <: AbstractRelation

The Fisher identity `γ = ν(2 − η)`, relating the susceptibility to the
correlation function that produces it.

Reference: [Fisher1964](@cite) (J. Math. Phys. **5**, 944).
"""
@relation :scaling Fisher(γ, ν, η) = γ - ν * (2 - η)

"""
    Josephson <: AbstractRelation

The Josephson (hyperscaling) identity `2 − α = d·ν`.  Valid below the
upper critical dimension; at and above it, mean-field exponents satisfy
it only at `d = d_upper` (e.g. `d = 4` for Ising).

Reference: [Josephson1967](@cite) (Proc. Phys. Soc. **92**, 269, "Inequality for
the specific heat: I. Derivation") — again an inequality first, `2 − α ≥ dν`.
"""
@relation :scaling Josephson(α, ν, d) = 2 - α - d * ν

"""
    DynamicalScaling <: AbstractRelation

The dynamical-exponent scaling of the gap with the correlation length on
approach to a quantum critical point,

`Δ ∼ ξ^{−z}`   ⟹   `d(ln Δ)/d(ln ξ) = −z`.

Supplied-derivative convention: `dlogΔ_dlogξ` is the caller-computed
log–log slope of the gap against the correlation length.  Reads the
dynamical critical exponent `z` off measured `(Δ, ξ)` pairs.

Holds only where a finite `z` exists; see [`ActivatedDynamicalScaling`](@ref)
for the infinite-randomness case, where none does.

Variables: `dlogΔ_dlogξ`, `z`.

Reference: [IgloiMonthus2005](@cite) Eq. (4.14), §4.1.3 — `t_r ∼ ξ^z`, the same
statement one inverse-time up.
"""
@relation :scaling DynamicalScaling(dlogΔ_dlogξ, z::DynamicalExponent) = dlogΔ_dlogξ + z

"""
    ActivatedDynamicalScaling <: AbstractRelation

Activated dynamic scaling at an infinite-randomness fixed point, where the gap
closes exponentially in a power of the length rather than as a power of it,

`ln(1/Δ) ∼ ξ^ψ`   ⟹   `d(ln[ln(1/Δ)])/d(ln ξ) = ψ`.

This is not [`DynamicalScaling`](@ref) with some other `z`: no finite `z`
describes such a point at all, since `−d(ln Δ)/d(ln ξ) = ψ·ln(1/Δ)` grows
without bound.

Supplied-derivative convention: `dloglogΔ_dlogξ` is the caller-computed slope of
`ln[ln(1/Δ)]` against `ln ξ`; `Δ < 1` is required for the inner log.

Variables: `dloglogΔ_dlogξ`, `ψ`.

Reference: [IgloiMonthus2005](@cite) Eq. (4.13), §4.1.3, which defines `ψ` by
`ln t_r ∼ ξ^ψ` and gives `ψ = 1/2` for the 1D random transverse-field Ising
chain; that value is Fisher's, [FisherDS1995](@cite).
"""
@relation :scaling ActivatedDynamicalScaling(dloglogΔ_dlogξ, ψ::ActivatedExponent) =
    dloglogΔ_dlogξ - ψ

"""
    ActivatedFiniteSizeScaling <: AbstractRelation

Finite-size scaling of a TYPICAL observable at an infinite-randomness fixed
point, on an open chain of length `L`:

`ln O_typ(L) ∼ −L^ψ`   ⟹   `d(ln[−ln O_typ])/d(ln L) = ψ`.

The counterpart of [`FiniteSizeGap`](@ref), and the contrast is the point: a
conformal critical point closes its finite-size gap as a POWER of `L`
(`2πvx/L`, periodic chain), an infinite-randomness one as a STRETCHED
EXPONENTIAL on an open one.  One sweep in `L` tells them apart.

This is not [`ActivatedDynamicalScaling`](@ref) restated.  That relation is
about the correlation length `ξ`, which diverges at criticality and so
constrains nothing at `δ = 0`; this one is about the system size, the only
scale left there.  The review states them as separate equations.

Supplied-derivative convention: `dloglogO_dlogL` is the caller-computed slope of
`ln[−ln O_typ]` against `ln L`; `O_typ < 1` is required for the inner log.

Variables: `dloglogO_dlogL`, `ψ`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.17), §A.3, which states it for any
infinite-disorder fixed point (`L^ψ ln m` as the scaling combination is §2.4).
The 1D instances are Eq. (4.12) for the gap and Eq. (4.6) for the surface
magnetization, both in §4.1 and both stated for free boundary conditions.
"""
@relation :scaling ActivatedFiniteSizeScaling(dloglogO_dlogL, ψ::ActivatedExponent) =
    dloglogO_dlogL - ψ

"""
    TypicalCorrelationLength <: AbstractRelation

At an infinite-randomness fixed point the typical correlation length is an
anomalous POWER of the average one,

`ξ_typ ∼ ξ^{1−ψ} ∼ |δ|^{−ν(1−ψ)}`   ⟹   `ν_typ = ν(1 − ψ)`.

So `C_typ(r)` still decays exponentially off criticality — it is the length that
is anomalous, not the functional form. `ν_typ < ν` for any `ψ > 0`: the typical
correlation length is the SMALLER of the two, and the two coincide only where
`ψ = 0`, i.e. where the fixed point is not infinite-randomness at all.

Reference: [IgloiMonthus2005](@cite) Eq. (9.4), §9.1.2 — derived in `d`
dimensions, and restated in the scaling-theory appendix just below Eq. (A.21).
Not a 1D statement.  The 1D values it is checked against are Table 1 (§4.1.2):
`ν = 2` is Eq. (4.9) and `ν_typ = 1` is Eq. (4.10), obtained separately.

Variables: `ν_typ`, `ν`, `ψ`.
"""
@relation :scaling TypicalCorrelationLength(ν_typ, ν, ψ::ActivatedExponent) =
    ν_typ - ν * (1 - ψ)

"""
    GriffithsExponentDivergence <: AbstractRelation

Why the Griffiths dynamical exponent varies continuously and why it stops
existing at the fixed point,

`z ∼ ξ^ψ ∼ |δ|^{−νψ}`   ⟹   `d(ln z)/d(ln|δ|) = −νψ`.

Off criticality `z(δ)` is finite and non-universal, so [`DynamicalScaling`](@ref)
holds with a `z` that drifts with distance from the transition; as `δ → 0` it
diverges, and [`ActivatedDynamicalScaling`](@ref) is what remains. The two
dynamic scalings are endpoints of one law, and `νψ` is the rate.

Supplied-derivative convention: `dlogz_dlogδ` is the caller-computed log–log
slope of `z` against `|δ|`.

Reference: [IgloiMonthus2005](@cite) Eq. (9.5), §9.1.2.  In 1D the Griffiths
exponent is known in closed form, `1/z = 2|δ|` (stated with Eq. (4.51),
§4.4.2), whose slope is `−1` — which is `−νψ` at that chain's `ν = 2`, `ψ = 1/2`.

Variables: `dlogz_dlogδ`, `ν`, `ψ`.
"""
@relation :scaling GriffithsExponentDivergence(dlogz_dlogδ, ν, ψ::ActivatedExponent) =
    dlogz_dlogδ + ν * ψ

"""
    GriffithsSusceptibility <: AbstractRelation

The Griffiths-phase susceptibility singularity, `χ(T) ∼ T^{−1+d/z}`, which is
what makes a continuously varying `z` measurable:

`d(ln χ)/d(ln T) = −1 + d/z`.

Divergent for `z > d`, finite for `z < d` — so the Griffiths phase has a line
inside it, at `z = d`, across which `χ` stops diverging while nothing else
happens. Written multiplied through by `z`, which keeps the residual affine in
every variable (so generic `solve` works) and finite at `z = 0`.

Reads two more of Appendix A's equations unchanged, because they carry the same
combination on other axes: the Griffiths gap DISTRIBUTION `P(ε) ∼ ε^{−1+d/z}`
(Eq. (A.29)), which is the microscopic origin of the rest, and the field-driven
`χ(H) ∼ H^{−1+d/z}` (Eq. (A.33)).

Reference: [IgloiMonthus2005](@cite) Eq. (4.56), §4.4.2, written there for the
1D chain, and Eq. (A.32), §A.4.1 in `d` dimensions.  §9.1.2 states the
`d`-dimensional form as the replacement `z → z/d` in Eqs. (4.51), (4.55) and
(4.56), which is the `d` carried here.

Variables: `dlogχ_dlogT`, `d`, `z`.
"""
@relation :scaling GriffithsSusceptibility(dlogχ_dlogT, d, z::DynamicalExponent) =
    (dlogχ_dlogT + 1) * z - d

"""
    GriffithsSpecificHeat <: AbstractRelation

The companion singularity in the same phase, `s(T) ∼ c_V(T) ∼ T^{d/z}`:

`d(ln c_V)/d(ln T) = d/z`.

Reads the same `z` as [`GriffithsSusceptibility`](@ref) off a different
observable, which is the point — a `z` fitted from one alone is a fit, and the
two together are a check. Multiplied through by `z`, as there.

Reference: [IgloiMonthus2005](@cite) Eq. (4.51), §4.4.2 (`s(T) ∼ c_V(T)`, so
this reads the entropy too), under the same `z → z/d` of §9.1.2.

Variables: `dlogc_dlogT`, `d`, `z`.
"""
@relation :scaling GriffithsSpecificHeat(dlogc_dlogT, d, z::DynamicalExponent) =
    dlogc_dlogT * z - d

"""
    ActivatedMomentGrowth <: AbstractRelation

The moment of a surviving cluster grows as a power of the LOG energy scale,
`μ ∼ |ln Ω|^φ`, with

`φ = (d − x_m)/ψ`.

This is where the fixed point's irrational exponents come from: it is not that
`φ` is separately measured, but that `ψ` converts a length dimension into a
log-energy one. For the 1D random transverse-field Ising chain
(`d = 1`, `x_m = (3−√5)/4`, `ψ = 1/2`) it returns the golden mean `(1+√5)/2`.

Written multiplied through by `ψ`, keeping the residual affine in every variable.

Reference: [IgloiMonthus2005](@cite) Eq. (A.21), §A.4, where `d − x_m` is
identified as the fractal dimension of the cluster.  The golden mean it returns
for the RTFIC is that review's Eq. (3.18), §3.5, reached by a different route.

Variables: `φ`, `d`, `x_m`, `ψ`.
"""
@relation :scaling ActivatedMomentGrowth(φ, d, x_m, ψ::ActivatedExponent) =
    φ * ψ - (d - x_m)

# ─── Appendix A: the four scaling types of a random system ───────────────
#
# [IgloiMonthus2005](@cite) Appendix A gives the same menu of observables for
# each kind of fixed point a random system can flow to: how the energy scale
# tracks the size, how autocorrelations decay, and what the low-temperature
# thermodynamics looks like.  The four answers are different in FORM, not just
# in exponent value, which is what makes measuring one of them a classification
# rather than a fit:
#
#   conventional random critical (A.2)  power of L,       power of t,   power of T
#   infinite disorder            (A.3)  exponential in L, power of ln t, power of ln T
#   Griffiths phase              (A.4)  power of L,       power of t,   power of T
#   large spin                   (A.5)  power of L,       -            , Curie
#
# The conventional and Griffiths columns agree in form and differ in which
# exponent combination appears, so the relations below are separate objects.
#
# What is NOT here: an equation whose exponent combination already appears
# above is the same relation on another axis, named in a docstring rather than
# duplicated.  Eqs. (A.26) and (A.33) restate observables against a small
# ordering field `H` instead of `T`; Eq. (A.29) is the Griffiths gap
# DISTRIBUTION, not a field statement at all, but carries the same `-1+d/z`.
# Only a NEW combination gets its own relation (Eq. (A.15)).

"""
    OrderParameterDimension <: AbstractRelation

The order-parameter exponent is its scaling dimension times the correlation
exponent, `β = ν·x_m`, which is what the finite-size and the thermodynamic
readings of one scaling form have to agree on: `m(δ) ∼ δ^β` as `L → ∞`, and
`m(δ=0, L) ∼ L^{−x_m}` at fixed `δ = 0`.

Holds verbatim at a SURFACE with `x_m` replaced by the surface dimension
`x_m^s`, giving `β_s = ν·x_m^s` (see [`SurfaceMagnetization`](@ref)). The two
readings are independent measurements, so this is a check and not a definition.

Variables: `β`, `ν`, `x_m`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.11) and the paragraph below it,
which states both the bulk and the surface form; the 1D instances are Table 1
(§4.1.2), where `β = νx_m` and `β_s = νx_m^s` are checked against each other.
"""
@relation :scaling OrderParameterDimension(β, ν, x_m::ScalingDimension) = β - ν * x_m

"""
    FixedPointDisorderStrength <: AbstractRelation

At a fixed point where the low-energy excitations are LOCALIZED, the disorder
strength stops flowing at a finite value tied to the dynamical exponent,

`D = z/d`.

This is the statement that separates the two random fixed points that are not
infinite-randomness ones. A conventional random critical point and a Griffiths
phase both reach it, from the same argument: the gap distribution `P(ε) ∼
ε^{−1+d/z}` is a fixed point of the rescaling only at this `D`. An
infinite-randomness point is exactly where it fails, `D` running away as the
scale grows, and that is what makes `z` infinite there and
[`ActivatedExponent`](@ref) the exponent that survives.

Variables: `D`, `z`, `d`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.10), §A.2 (conventional random
critical scaling) and restated below Eq. (A.29), §A.4.1, for the Griffiths
phase. The contrast is Eq. (A.18), §A.3.
"""
@relation :scaling FixedPointDisorderStrength(
    D::DisorderStrength, z::DynamicalExponent, d
) = D * d - z

"""
    ConventionalFiniteSizeEnergy <: AbstractRelation

The finite-size energy scale of a random system whose excitations are localized,

`Ω ∼ L^{−z}`   ⟹   `d(ln Ω)/d(ln L) = −z`.

The size-space counterpart of [`DynamicalScaling`](@ref), and the form
[`ActivatedFiniteSizeScaling`](@ref) replaces. Unlike [`FiniteSizeGap`](@ref)
it carries no amplitude, so it applies where there is no conformal field theory
to supply one: any `z`, periodic or open.

Variables: `dlogΩ_dlogL`, `z`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.8), §A.2, restated as Eq. (A.30),
§A.4.1 for the Griffiths phase and as Eq. (A.38), §A.5 for the large-spin
phase. One form, three fixed points.
"""
@relation :scaling ConventionalFiniteSizeEnergy(dlogΩ_dlogL, z::DynamicalExponent) =
    dlogΩ_dlogL + z

"""
    CriticalAutocorrelation <: AbstractRelation

Autocorrelation decay at a conventional random critical point,

`G(t) ∼ t^{−2x_m/z}`   ⟹   `d(ln G)/d(ln t)·z + 2x_m = 0`.

Both exponents enter, so a measured decay reads neither alone: pair it with
[`ConventionalFiniteSizeEnergy`](@ref) for `z`, or with
[`OrderParameterDimension`](@ref) for `x_m`. Multiplied through by `z` to stay
affine in each variable.

Variables: `dlogG_dlogt`, `x_m`, `z`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.13), §A.2. Contrast
[`ActivatedAutocorrelation`](@ref), where `t` is replaced by `ln t`, and
[`GriffithsAutocorrelation`](@ref), where `2x_m` is replaced by `d`.
"""
@relation :scaling CriticalAutocorrelation(
    dlogG_dlogt, x_m::ScalingDimension, z::DynamicalExponent
) = dlogG_dlogt * z + 2 * x_m

"""
    CriticalQuantumSusceptibility <: AbstractRelation

Low-temperature susceptibility at a conventional random QUANTUM critical point,

`χ(T) ∼ T^{−γ/νz}`   ⟹   `d(ln χ)/d(ln T)·νz + γ = 0`,

from `T` setting an energy scale and hence a thermal length `L_T ∼ T^{−1/z}`.
The `z` is what makes this a quantum statement: at `z = 1` it is the classical
`χ ∼ |t|^{−γ}` read along the temperature axis.

Variables: `dlogχ_dlogT`, `γ`, `ν`, `z`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.14), §A.2. The Griffiths-phase
counterpart is [`GriffithsSusceptibility`](@ref) and the infinite-randomness
one is [`ActivatedSusceptibility`](@ref); all three are different forms, not
different values.
"""
@relation :scaling CriticalQuantumSusceptibility(dlogχ_dlogT, γ, ν, z::DynamicalExponent) =
    dlogχ_dlogT * ν * z + γ

"""
    CriticalQuantumSpecificHeat <: AbstractRelation

Low-temperature specific heat at a conventional random quantum critical point,

`c_V(T) ∼ T^{−α/νz}`   ⟹   `d(ln c_V)/d(ln T)·νz + α = 0`.

Same thermal length as [`CriticalQuantumSusceptibility`](@ref), so the two
together over-determine `z` from thermodynamics alone.

Variables: `dlogc_dlogT`, `α`, `ν`, `z`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.14), §A.2.
"""
@relation :scaling CriticalQuantumSpecificHeat(dlogc_dlogT, α, ν, z::DynamicalExponent) =
    dlogc_dlogT * ν * z + α

"""
    ActivatedAutocorrelation <: AbstractRelation

Autocorrelation decay at an infinite-randomness fixed point is ULTRA-SLOW,
logarithmic in time rather than a power of it,

`G(t) ∼ (ln t)^{−x_m/ψ}`   ⟹   `d(ln G)/d(ln ln t)·ψ + x_m = 0`.

Not a small exponent in [`CriticalAutocorrelation`](@ref): no power of `t`
describes this at all. The reason is that disorder is strictly correlated along
the time direction, so a region that is locally ordered at time 0 stays ordered,
and the density of such rare regions is what decays.

Supplied-derivative convention: the slope of `ln G` against `ln ln t`, so `t > e`
is required.

Variables: `dlogG_dloglnt`, `x_m`, `ψ`.

The field-driven susceptibility carries the same `x_m/ψ` on the `ln H` axis,
`χ(H)·H ∼ (ln H)^{−x_m/ψ}` (Eq. (A.26)), so this relation reads that too.

Reference: [IgloiMonthus2005](@cite) Eq. (A.23), §A.3, derived from the scaling
form Eq. (A.22).
"""
@relation :scaling ActivatedAutocorrelation(
    dlogG_dloglnt, x_m::ScalingDimension, ψ::ActivatedExponent
) = dlogG_dloglnt * ψ + x_m

"""
    ActivatedSusceptibility <: AbstractRelation

Low-temperature susceptibility at an infinite-randomness fixed point: a Curie
law times a POWER OF THE LOG of the temperature,

`χ(T) ∼ (ln T)^{(d−2x_m)/ψ} / T`   ⟹   `d(ln[χT])/d(ln|ln T|)·ψ − (d − 2x_m) = 0`.

The rare regions each contribute a Curie term, so the `1/T` is not a critical
singularity at all and has to be divided out before the exponent is read. Doing
that is what the supplied-derivative convention says: the slope of `ln(χT)`
against `ln|ln T|`.

Variables: `dlogχT_dloglnT`, `d`, `x_m`, `ψ`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.25), §A.3.
"""
@relation :scaling ActivatedSusceptibility(
    dlogχT_dloglnT, d, x_m::ScalingDimension, ψ::ActivatedExponent
) = dlogχT_dloglnT * ψ - (d - 2 * x_m)

"""
    ActivatedSpecificHeat <: AbstractRelation

Low-temperature specific heat at an infinite-randomness fixed point,

`c_V(T) ∼ (ln T)^{−d/ψ}`   ⟹   `d(ln c_V)/d(ln|ln T|)·ψ + d = 0`,

from the rare low-energy excitations being a distance `L_T ∼ (ln T)^{1/ψ}` apart,
so `c_V ∼ L_T^{−d}`. Note what is absent: no exponent of the ordered phase
enters, only `d` and `ψ`, so this reads `ψ` off thermodynamics with nothing
else fitted.

Unchanged on the field axis, `c_V(H) ∼ (ln H)^{−d/ψ}` (Eq. (A.26)).

Variables: `dlogc_dloglnT`, `d`, `ψ`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.25), §A.3.
"""
@relation :scaling ActivatedSpecificHeat(dlogc_dloglnT, d, ψ::ActivatedExponent) =
    dlogc_dloglnT * ψ + d

"""
    GriffithsAutocorrelation <: AbstractRelation

Autocorrelation decay in a Griffiths phase,

`G(t) ∼ t^{−d/z}`   ⟹   `d(ln G)/d(ln t)·z + d = 0`.

A power law with a CONTINUOUSLY VARYING exponent: `z = z(δ)` depends on the
distance from criticality, so unlike every other relation here the exponent is
not universal, only the form is. Averaging `exp(−t/t_r)` over the algebraic tail
`p(t_r) ∼ t_r^{−d/z−1}` of rare-region relaxation times is where it comes from,
which is also why the decay is set by `d` rather than by `x_m`.

Reads the same `z` as [`GriffithsSusceptibility`](@ref) and
[`GriffithsSpecificHeat`](@ref), off dynamics rather than thermodynamics.

Variables: `dlogG_dlogt`, `d`, `z`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.27), §A.4.1. In the ORDERED
Griffiths phase of a chain the same argument gives `2/z` instead of `d/z`,
because isolating a domain costs two weak bonds rather than a surface.
"""
@relation :scaling GriffithsAutocorrelation(dlogG_dlogt, d, z::DynamicalExponent) =
    dlogG_dlogt * z + d

"""
    LargeSpinMoment <: AbstractRelation

The large-spin fixed point of a random chain with mixed ferromagnetic and
antiferromagnetic couplings, where renormalization grows an effective spin
rather than freezing singlets: `S_eff ∼ L^{dζ}` with `Ω ∼ L^{−z}` gives

`S_eff ∼ Ω^{−κ}`,  `κ = dζ/z`.

The fourth kind of random fixed point, and the only one whose order parameter
GROWS under renormalization. `ζ = 1/2` follows from a random-walk argument on
the signs of the couplings.

Variables: `κ`, `d`, `ζ`, `z`.

Reference: [IgloiMonthus2005](@cite) Eqs. (A.37) and (A.38), §A.5; the 1D
instance is Eqs. (8.6)-(8.8), §8.2, where `ζ = 1/2` and `κ = 0.22(1)` is
measured, giving `z = 1/(2κ)`.
"""
@relation :scaling LargeSpinMoment(κ, d, ζ, z::DynamicalExponent) = κ * z - d * ζ

"""
    ConventionalFieldSusceptibility <: AbstractRelation

Susceptibility of a conventional random quantum critical point against a small
ORDERING FIELD rather than against temperature,

`χ(H) ∼ H^{−γ/[ν(d+z−x_m)]}`.

The denominator is the whole difference from
[`CriticalQuantumSusceptibility`](@ref): a field couples to the order parameter
over a correlation volume, so the length it sets is `L_H ∼ H^{−1/(d+z−x_m)}`
rather than the thermal `L_T ∼ T^{−1/z}`.  The two are independent readings of
`z`, and `x_m` enters only the field one.

Variables: `dlogχ_dlogH`, `γ`, `ν`, `d`, `z`, `x_m`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.15), §A.2.
"""
@relation :scaling ConventionalFieldSusceptibility(
    dlogχ_dlogH, γ, ν, d, z::DynamicalExponent, x_m::ScalingDimension
) = dlogχ_dlogH * ν * (d + z - x_m) + γ

"""
    ConventionalFieldSpecificHeat <: AbstractRelation

Specific heat of a conventional random quantum critical point against a small
ordering field,

`c_V(H) ∼ H^{−α/[ν(d+z−x_m)]}`.

Same field length as [`ConventionalFieldSusceptibility`](@ref), so the pair
over-determines `d + z − x_m` from field sweeps alone, exactly as
[`CriticalQuantumSusceptibility`](@ref) and
[`CriticalQuantumSpecificHeat`](@ref) over-determine `νz` from temperature ones.

Variables: `dlogc_dlogH`, `α`, `ν`, `d`, `z`, `x_m`.

Reference: [IgloiMonthus2005](@cite) Eq. (A.15), §A.2.
"""
@relation :scaling ConventionalFieldSpecificHeat(
    dlogc_dlogH, α, ν, d, z::DynamicalExponent, x_m::ScalingDimension
) = dlogc_dlogH * ν * (d + z - x_m) + α

"""
    OrderedGriffithsEnergyScale <: AbstractRelation

The ORDERED Griffiths phase above one dimension has an energy scale that is
neither a power of the size nor a stretched exponential in it, but a power of
the LOG of it,

`|ln Ω| ∼ (ln L)^{1/d}`   ⟹   `d(ln|ln Ω|)/d(ln ln L)·d = 1`.

A fifth functional form, and the reason the ordered and the disordered
Griffiths phases are not mirror images: isolating an ordered cluster of `l^d`
sites needs a moat of width `∼ l^d`, so the cost goes as `l^{d²}` and the
rare-region statistics change shape.  At `d = 1` the slope is 1 and the law
collapses back to [`ConventionalFiniteSizeEnergy`](@ref), a plain power of `L`.

Variables: `dloglogΩ_dloglogL`, `d`.

The same `|ln Ω| ∼ (ln L)^x` form turns up elsewhere with an exponent that is
not `1/d` (§9.3 quotes `x = 2` and `x = 3/2` as conjectures for a 2D Dirac
problem), where it reads as a divergent `z` or a vanishing `ψ`.  This relation
is the case where `x` is derived rather than fitted.

Reference: [IgloiMonthus2005](@cite) Eq. (A.35), §A.4.2, with the companion
autocorrelation `G(t) ∼ exp(−A|ln t|^d)` of Eq. (A.34) (= Eq. (9.7), §9.1.3).
"""
@relation :scaling OrderedGriffithsEnergyScale(dloglogΩ_dloglogL, d) =
    dloglogΩ_dloglogL * d - 1

"""
    WeinribHalperinExponent <: AbstractRelation

When spatially correlated disorder IS relevant, the correlation-length exponent
it flows to is fixed by the correlation decay alone,

`ν = 2/ρ`   for   `G_d(r) ∼ r^{−ρ}` with `ρ < 2/ν_unc`,

where `ν_unc` is the correlation-length exponent of the same model with
UNCORRELATED disorder, not the clean `ν₀` that [`HarrisCriterion`](@ref) takes.

The companion of [`WeinribHalperinCriterion`](@ref), which decides whether that
`ρ` matters at all: the criterion's marginal line `ρ = 2/ν` and this relation
are the same equation, so the new exponent is exactly where the perturbation
stops being relevant.  Nothing about the clean model survives except through the
threshold, which is why correlated disorder makes a new universality class
rather than shifting the old one.

Not a 1D statement.  The review derives it for the random transverse-field
chain and then states it holds in higher dimensions too, on the general argument
of Weinrib and Halperin ([WeinribHalperin1983](@cite)).

Variables: `ν_dis`, `ρ`.

Reference: [IgloiMonthus2005](@cite) Eq. (10.3), §10.1 for the 1D case and the
paragraph below it for the general one; the relevance threshold is Eq. (10.2).
"""
@relation :scaling WeinribHalperinExponent(ν_dis, ρ::DisorderCorrelationExponent) =
    ν_dis * ρ - 2

"""
    exponents_consistent(nt::NamedTuple; d, atol=0) -> Bool

Gate-check a `CriticalExponents`-style NamedTuple `(α, β, γ, δ, ν, η)`
against all scaling relations at spatial dimension `d`.  A thin wrapper
over the generic registry sweep ([`check_all`](@ref) with
`domain = :scaling`); kept as the domain-specific entry point atlases
gate their exponent tables with.
"""
function exponents_consistent(nt::NamedTuple; d, atol=0)
    return check_all((; nt..., d=d); atol=atol, domain=:scaling)
end
export exponents_consistent

"""
    exponents_consistent(s::ScalingDimensions; atol=0) -> Bool
    exponents_consistent(s::InfiniteRandomness; atol=0) -> Bool

The same gate run on a fixed point that carries its own `d`, so the dimension
cannot be restated wrongly at the call site.  Equivalent to
`exponents_consistent(critical_exponents(s); d = s.d)`.

This is the door to prefer.  Passing a bare NamedTuple leaves `d` to the caller,
and a sweep silently SKIPS every relation whose variables are not all present
(see [`applicable_relations`](@ref)), so forgetting `d` does not fail; it just
stops checking the twelve relations that need it.
"""
function exponents_consistent(s::ScalingDimensions; atol=0)
    return exponents_consistent(critical_exponents(s); d=s.d, atol=atol)
end
function exponents_consistent(s::InfiniteRandomness; atol=0)
    return exponents_consistent(critical_exponents(s); d=s.d, atol=atol)
end

"""
    exponent_residuals(s::ScalingDimensions) -> NamedTuple
    exponent_residuals(s::InfiniteRandomness) -> NamedTuple

Per-relation residuals of a fixed point that carries its own `d`, the
diagnostic companion of [`exponents_consistent`](@ref).
"""
exponent_residuals(s::ScalingDimensions) = exponent_residuals(critical_exponents(s); d=s.d)
exponent_residuals(s::InfiniteRandomness) = exponent_residuals(critical_exponents(s); d=s.d)

"""
    exponent_residuals(nt::NamedTuple; d) -> NamedTuple

Per-relation residuals of the scaling laws for the exponent set
`nt = (α, β, γ, δ, ν, η)` at dimension `d` — the diagnostic companion of
[`exponents_consistent`](@ref), keyed by lowercase relation name.
"""
function exponent_residuals(nt::NamedTuple; d)
    report = relation_report((; nt..., d=d); domain=:scaling)
    names = Tuple(Symbol(lowercase(String(nameof(typeof(row.relation))))) for row in report)
    return NamedTuple{names}(Tuple(row.residual for row in report))
end
export exponent_residuals
