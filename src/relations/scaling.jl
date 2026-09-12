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

Reference: [IgloiMonthus2005](@cite) Eq. (4.12) for the gap and Eq. (4.6) for
the surface magnetization, both in §4.1 and both stated for free boundary
conditions; `L^ψ ln m` as the scaling combination of any infinite-disorder
fixed point is §2.4.
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

Reference: [IgloiMonthus2005](@cite) Eq. (4.56), §4.4.2, written there for the
1D chain.  §9.1.2 states the `d`-dimensional form as the replacement `z → z/d`
in Eqs. (4.51), (4.55) and (4.56), which is the `d` carried here.

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
