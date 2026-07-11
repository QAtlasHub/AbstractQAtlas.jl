# relations/scaling.jl — critical-exponent scaling laws as relations.
#
# The four classical identities among the equilibrium exponents
# {α, β, γ, δ, ν, η} (+ dimension d for hyperscaling).  Exact exponent
# sets (2D Ising rationals, mean-field at the upper critical dimension)
# satisfy them with residual ≡ 0 in exact arithmetic; numerical sets
# (3D Ising bootstrap) satisfy them within quoted uncertainties.
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
struct Rushbrooke <: AbstractRelation end
export Rushbrooke

residual(::Rushbrooke; α, β, γ) = α + 2β + γ - 2
solve(::Rushbrooke, ::Val{:α}; β, γ) = 2 - 2β - γ
solve(::Rushbrooke, ::Val{:β}; α, γ) = (2 - α - γ) / 2
solve(::Rushbrooke, ::Val{:γ}; α, β) = 2 - α - 2β

"""
    Widom <: AbstractRelation

The Widom identity `γ = β(δ − 1)`.
"""
struct Widom <: AbstractRelation end
export Widom

residual(::Widom; β, γ, δ) = γ - β * (δ - 1)
solve(::Widom, ::Val{:β}; γ, δ) = γ / (δ - 1)
solve(::Widom, ::Val{:γ}; β, δ) = β * (δ - 1)
solve(::Widom, ::Val{:δ}; β, γ) = γ / β + 1

"""
    Fisher <: AbstractRelation

The Fisher identity `γ = ν(2 − η)`.
"""
struct Fisher <: AbstractRelation end
export Fisher

residual(::Fisher; γ, ν, η) = γ - ν * (2 - η)
solve(::Fisher, ::Val{:γ}; ν, η) = ν * (2 - η)
solve(::Fisher, ::Val{:ν}; γ, η) = γ / (2 - η)
solve(::Fisher, ::Val{:η}; γ, ν) = 2 - γ / ν

"""
    Josephson <: AbstractRelation

The Josephson (hyperscaling) identity `2 − α = d·ν`.  Valid below the
upper critical dimension; at and above it, mean-field exponents satisfy
it only at `d = d_upper` (e.g. `d = 4` for Ising).
"""
struct Josephson <: AbstractRelation end
export Josephson

residual(::Josephson; α, ν, d) = 2 - α - d * ν
solve(::Josephson, ::Val{:α}; ν, d) = 2 - d * ν
solve(::Josephson, ::Val{:ν}; α, d) = (2 - α) / d
solve(::Josephson, ::Val{:d}; α, ν) = (2 - α) / ν

"""
    exponents_consistent(nt::NamedTuple; d, atol=0) -> Bool

Gate-check a `CriticalExponents`-style NamedTuple `(α, β, γ, δ, ν, η)`
against all four scaling relations at spatial dimension `d`.  With the
default `atol = 0` this demands exact satisfaction (appropriate for
`Rational` exponent sets); pass a finite `atol` for numerical sets.

Intended use: an atlas registers an exponent table and gates it with
this function, instead of re-deriving individual relations in comments.

See also [`exponent_residuals`](@ref) for per-relation diagnostics.
"""
function exponents_consistent(nt::NamedTuple; d, atol=0)
    return all(abs(r) <= atol for r in values(exponent_residuals(nt; d=d)))
end
export exponents_consistent

"""
    exponent_residuals(nt::NamedTuple; d) -> NamedTuple

Per-relation residuals of the four scaling laws for the exponent set
`nt = (α, β, γ, δ, ν, η)` at dimension `d` — the diagnostic companion of
[`exponents_consistent`](@ref).
"""
function exponent_residuals(nt::NamedTuple; d)
    return (
        rushbrooke=residual(Rushbrooke(); α=nt.α, β=nt.β, γ=nt.γ),
        widom=residual(Widom(); β=nt.β, γ=nt.γ, δ=nt.δ),
        fisher=residual(Fisher(); γ=nt.γ, ν=nt.ν, η=nt.η),
        josephson=residual(Josephson(); α=nt.α, ν=nt.ν, d=d),
    )
end
export exponent_residuals
