# relations/scaling.jl — critical-exponent scaling laws.
#
# Each law is one @relation declaration; residual/check/solve for every
# variable follow mechanically (all four laws are affine in every
# variable).  Exact exponent sets (2D Ising rationals, mean-field at the
# upper critical dimension) satisfy them with residual ≡ 0 in exact
# arithmetic; numerical sets (3D Ising bootstrap) within quoted errors.
#
# References (textbook standard): Rushbrooke, [Rushbrooke1963](@cite); Widom,
# [Widom1965](@cite); Fisher, [Fisher1964](@cite); Josephson,
# [Josephson1967](@cite).  The infinite-randomness block below is
# Iglói–Monthus, [IgloiMonthus2005](@cite), whose 1D values are Fisher,
# [FisherDS1995](@cite).

"""
    Rushbrooke <: AbstractRelation

The Rushbrooke identity `α + 2β + γ = 2`.

```julia
residual(Rushbrooke(); α=0//1, β=1//8, γ=7//4)   # == 0//1 (2D Ising, exact)
solve(Rushbrooke(), Val(:γ); α=0//1, β=1//8)     # == 7//4
```
"""
@relation :scaling Rushbrooke(α, β, γ) = α + 2β + γ - 2

"""
    Widom <: AbstractRelation

The Widom identity `γ = β(δ − 1)`.
"""
@relation :scaling Widom(β, γ, δ) = γ - β * (δ - 1)

"""
    Fisher <: AbstractRelation

The Fisher identity `γ = ν(2 − η)`.
"""
@relation :scaling Fisher(γ, ν, η) = γ - ν * (2 - η)

"""
    Josephson <: AbstractRelation

The Josephson (hyperscaling) identity `2 − α = d·ν`.  Valid below the
upper critical dimension; at and above it, mean-field exponents satisfy
it only at `d = d_upper` (e.g. `d = 4` for Ising).
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
"""
@relation :scaling ActivatedDynamicalScaling(dloglogΔ_dlogξ, ψ::ActivatedExponent) =
    dloglogΔ_dlogξ - ψ

"""
    TypicalCorrelationLength <: AbstractRelation

At an infinite-randomness fixed point the typical correlation length is an
anomalous POWER of the average one,

`ξ_typ ∼ ξ^{1−ψ} ∼ |δ|^{−ν(1−ψ)}`   ⟹   `ν_typ = ν(1 − ψ)`.

So `C_typ(r)` still decays exponentially off criticality — it is the length that
is anomalous, not the functional form. `ν_typ < ν` for any `ψ > 0`: the typical
correlation length is the SMALLER of the two, and the two coincide only where
`ψ = 0`, i.e. where the fixed point is not infinite-randomness at all.

Iglói–Monthus, [IgloiMonthus2005](@cite) Eq. (C_typ), derived in `d` dimensions
and restated in the general scaling section — this is not a 1D statement.

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

Iglói–Monthus, [IgloiMonthus2005](@cite) Eq. (z_delta).

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

Iglói–Monthus, [IgloiMonthus2005](@cite) Eq. (chi_T), with the stated
`d`-dimensional replacement `z → z/d`.

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

Iglói–Monthus, [IgloiMonthus2005](@cite) Eq. (entropy_d), with `z → z/d`.

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

Iglói–Monthus, [IgloiMonthus2005](@cite) Eq. (scaling_mu2).

Variables: `φ`, `d`, `x_m`, `ψ`.
"""
@relation :scaling ActivatedMomentGrowth(φ, d, x_m, ψ::ActivatedExponent) =
    φ * ψ - (d - x_m)

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
