# relations/scaling.jl — critical-exponent scaling laws.
#
# Each law is one @relation declaration; residual/check/solve for every
# variable follow mechanically (all four laws are affine in every
# variable).  Exact exponent sets (2D Ising rationals, mean-field at the
# upper critical dimension) satisfy them with residual ≡ 0 in exact
# arithmetic; numerical sets (3D Ising bootstrap) within quoted errors.
#
# References (textbook standard): Rushbrooke, J. Chem. Phys. 39, 842
# (1963); Widom, J. Chem. Phys. 43, 3898 (1965); Fisher, J. Math. Phys.
# 5, 944 (1964); Josephson, Proc. Phys. Soc. 92, 269 (1967).

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
