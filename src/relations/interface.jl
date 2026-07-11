# relations/interface.jl — the uniform three-verb contract for relations.
#
# A relation is an *identity among observables or exponents* — a generic,
# model-independent statement (a scaling law, a fluctuation–dissipation
# identity, a topological quantization).  Expressing each one ONCE, as a
# tested first-class object, is the point of this package: downstream
# atlases stop re-deriving them ad hoc in comments and per-model tests.

"""
    AbstractRelation

Abstract parent type for physics relations.  Every concrete relation is
a singleton struct implementing (a subset of) the three verbs:

- [`residual`](@ref)`(rel; vars...)` — signed violation; `0` ⇔ satisfied.
- [`check`](@ref)`(rel; atol=0, vars...)` — `|residual| ≤ atol`.
- [`solve`](@ref)`(rel, Val(:x); vars...)` — the value of `x` implied by
  the remaining variables.

**Exact-arithmetic contract**: `residual` and `solve` must not promote
their inputs — `Rational` in ⇒ `Rational` out, so exactly-known values
(e.g. the 2D Ising exponents) satisfy their relations *exactly*
(`residual == 0//1`), not merely to floating-point tolerance.
"""
abstract type AbstractRelation end
export AbstractRelation

"""
    residual(rel::AbstractRelation; vars...) -> Number

Signed violation of the relation at the given variable values; zero if
and only if the relation is satisfied.  Preserves the input number types
(see the exact-arithmetic contract on [`AbstractRelation`](@ref)).
"""
function residual end
export residual

"""
    check(rel::AbstractRelation; atol=0, vars...) -> Bool

`abs(residual(rel; vars...)) ≤ atol`.  With the default `atol = 0` this
is an *exact* test — appropriate for `Rational` inputs; pass an explicit
`atol` for floating-point / error-bar data.
"""
function check(rel::AbstractRelation; atol=0, kwargs...)
    return abs(residual(rel; kwargs...)) <= atol
end
export check

"""
    solve(rel::AbstractRelation, ::Val{x}; vars...) -> Number

The value of variable `x` implied by the relation and the remaining
variables, e.g. `solve(Widom(), Val(:γ); β=1//8, δ=15//1) == 7//4`.
Implemented per-relation for each solvable variable; preserves input
number types.
"""
function solve end
export solve
