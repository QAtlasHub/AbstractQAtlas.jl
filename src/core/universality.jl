# core/universality.jl — universality-class dispatch machinery.
#
# Ported from QAtlas `src/core/universality.jl`: the parametric dispatch
# tag and the quantity tags it answers.  Per-class DATA (critical
# exponents, central charges) and universal-BEHAVIOUR implementations
# stay in the implementing atlas — this package only owns the vocabulary
# plus the model-independent *relations between* exponents
# (see `relations/scaling.jl`).

"""
    Universality{C}

Parametric dispatch tag for universality classes. `C` is a `Symbol`
identifying the class (`:Ising`, `:XY`, `:Heisenberg`, `:Potts3`,
`:Potts4`, `:Percolation`, `:KPZ`, etc.).

Use with [`CriticalExponents`](@ref) (equilibrium) or
[`GrowthExponents`](@ref) (KPZ-type) and a `d` keyword to select the
spatial dimension:

```julia
fetch(Universality(:Ising), CriticalExponents(); d=2)   # exact Rational
fetch(Universality(:Ising), CriticalExponents(); d=3)   # numerical + _err
```
"""
struct Universality{C} <: AbstractQAtlasModel end
Universality(name::Symbol) = Universality{name}()
export Universality

"""
    CriticalExponents() <: AbstractQuantity

Standard set of equilibrium critical exponents
{α, β, γ, δ, ν, η} of a universality class. Returns a `NamedTuple`.

For exact values: fields are `Rational{Int}`.
For numerical estimates: fields are `Float64` with corresponding
`_err` fields (e.g., `β_err`) giving the uncertainty.

The scaling relations these exponents must satisfy are first-class
objects in this package — see [`Rushbrooke`](@ref), [`Widom`](@ref),
[`Fisher`](@ref), [`Josephson`](@ref) and the convenience gate
[`exponents_consistent`](@ref).
"""
struct CriticalExponents <: AbstractQuantity end
export CriticalExponents

"""
    GrowthExponents() <: AbstractQuantity

KPZ-type growth / roughness / dynamic exponents.  Returns
`(β_growth, α_rough, z)` instead of the equilibrium set.
"""
struct GrowthExponents <: AbstractQuantity end
export GrowthExponents
