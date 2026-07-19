# structure/potentials.jl — the thermodynamic-potential structure, and the Maxwell
# relations DERIVED from it.
#
# Each potential Φ(x, y) is a function of two natural variables; its differential is
#   dΦ = sₓ·cₓ dx + s_y·c_y dy,
# where the conjugate coefficient c (a state variable) is paired with the natural
# variable it multiplies and carried with its sign s = ±1.  The EQUALITY OF THE MIXED
# SECOND PARTIALS — ∂²Φ/∂x∂y = ∂²Φ/∂y∂x — is ONE structural fact; the four Maxwell
# relations are its four instances (U, F, H, G), not four independent axioms.  This
# layer declares the four differentials ONCE and DERIVES each Maxwell relation from
# ∂(sₓcₓ)/∂y = ∂(s_yc_y)/∂x — nothing hand-restated (the same move
# structure/criticality.jl makes for the scaling laws).

"""
    PotentialTerm(variable, conjugate, sign)

One term of a thermodynamic potential's differential:
`dΦ ⊃ sign · conjugate · d(variable)`.  The state variable `conjugate` is conjugate to
the natural `variable`, carried with its `sign` (`+1`/`−1`).
"""
struct PotentialTerm
    variable::Symbol       # the natural variable, e.g. :T
    conjugate::Symbol      # its conjugate state variable, e.g. :S
    sign::Int              # ±1, the sign of the term in dΦ
end

"""
    ThermodynamicPotential(name, x::PotentialTerm, y::PotentialTerm)

A thermodynamic potential `Φ` of the two natural variables `x.variable`, `y.variable`,
with differential
`dΦ = x.sign·x.conjugate·d(x.variable) + y.sign·y.conjugate·d(y.variable)`.  The four
standard potentials are [`thermodynamic_potentials`](@ref); the Maxwell relation falls
out via [`maxwell_relation`](@ref).
"""
struct ThermodynamicPotential
    name::Symbol
    x::PotentialTerm
    y::PotentialTerm
end

function Base.show(io::IO, p::ThermodynamicPotential)
    trm(c) = string(c.sign < 0 ? "−" : "+", c.conjugate, " d", c.variable)
    print(
        io,
        "ThermodynamicPotential(",
        p.name,
        "(",
        p.x.variable,
        ",",
        p.y.variable,
        "): d",
        p.name,
        " = ",
        trm(p.x),
        " ",
        trm(p.y),
        ")",
    )
    return nothing
end

"""
    thermodynamic_potentials() -> NTuple{4,ThermodynamicPotential}

The four standard thermodynamic potentials with their differentials — the single
structural source the Maxwell relations are derived from:

| Φ        | dΦ             |
|:---------|:---------------|
| `U(S,V)` | `+T dS − p dV` |
| `F(T,V)` | `−S dT − p dV` |
| `H(S,p)` | `+T dS + V dp` |
| `G(T,p)` | `−S dT + V dp` |
"""
function thermodynamic_potentials()
    return (
        ThermodynamicPotential(:U, PotentialTerm(:S, :T, +1), PotentialTerm(:V, :p, -1)),  # dU = +T dS − p dV
        ThermodynamicPotential(:F, PotentialTerm(:T, :S, -1), PotentialTerm(:V, :p, -1)),  # dF = −S dT − p dV
        ThermodynamicPotential(:H, PotentialTerm(:S, :T, +1), PotentialTerm(:p, :V, +1)),  # dH = +T dS + V dp
        ThermodynamicPotential(:G, PotentialTerm(:T, :S, -1), PotentialTerm(:p, :V, +1)),  # dG = −S dT + V dp
    )
end

export PotentialTerm, ThermodynamicPotential, thermodynamic_potentials

"""
    MaxwellRelation(potential, lhs, rhs, coeff)

The Maxwell relation DERIVED from a potential: `∂cₓ/∂y = coeff · ∂c_y/∂x`, where
`lhs = (cₓ, y)` names the derivative `∂cₓ/∂y`, `rhs = (c_y, x)` names `∂c_y/∂x`, and
`coeff = sₓ·s_y ∈ {+1,−1}`.  Its residual (zero ⇔ satisfied) is `∂cₓ/∂y − coeff·∂c_y/∂x`.
Built by [`maxwell_relation`](@ref).
"""
struct MaxwellRelation
    potential::Symbol
    lhs::Tuple{Symbol,Symbol}   # (conjugate, variable) of ∂cₓ/∂y
    rhs::Tuple{Symbol,Symbol}   # (conjugate, variable) of ∂c_y/∂x
    coeff::Int                  # ±1
end

function Base.show(io::IO, m::MaxwellRelation)
    eq = m.coeff < 0 ? " = −" : " = "
    print(
        io,
        "MaxwellRelation(",
        m.potential,
        ": ∂",
        m.lhs[1],
        "/∂",
        m.lhs[2],
        eq,
        "∂",
        m.rhs[1],
        "/∂",
        m.rhs[2],
        ")",
    )
    return nothing
end

"""
    maxwell_relation(p::ThermodynamicPotential) -> MaxwellRelation

Derive the Maxwell relation from a potential's differential.  Because `∂Φ/∂x = sₓ·cₓ`
and `∂Φ/∂y = s_y·c_y`, the commuting mixed partials `∂²Φ/∂x∂y = ∂²Φ/∂y∂x` give
`sₓ·∂cₓ/∂y = s_y·∂c_y/∂x`, i.e. `∂cₓ/∂y = (sₓ·s_y)·∂c_y/∂x` — the four Maxwell relations
are the four instances of this ONE identity.

```julia
maxwell_relation(thermodynamic_potentials()[2])   # F: ∂S/∂V = ∂p/∂T
```
"""
function maxwell_relation(p::ThermodynamicPotential)
    return MaxwellRelation(
        p.name,
        (p.x.conjugate, p.y.variable),   # ∂cₓ/∂y
        (p.y.conjugate, p.x.variable),   # ∂c_y/∂x
        p.x.sign * p.y.sign,
    )
end

"""
    maxwell_residual(m::MaxwellRelation; derivs) -> Number
    maxwell_residual(p::ThermodynamicPotential; derivs) -> Number

The residual of the structure-derived Maxwell relation, `∂cₓ/∂y − coeff·∂c_y/∂x`, with
the two first-derivative values read from `derivs` (an `AbstractDict` keyed by the
`(conjugate, variable)` tuples — the same keys `m.lhs` / `m.rhs` carry).  Zero ⇔ the
mixed partials of Φ commute.  Preserves the input number type (`Rational` in ⇒
`Rational` out).

```julia
m = maxwell_relation(thermodynamic_potentials()[2])              # F: ∂S/∂V = ∂p/∂T
maxwell_residual(m; derivs = Dict((:S, :V) => 2 // 1, (:p, :T) => 2 // 1))   # 0//1
```
"""
function maxwell_residual(m::MaxwellRelation; derivs)
    haskey(derivs, m.lhs) ||
        error("maxwell_residual: missing derivative ∂$(m.lhs[1])/∂$(m.lhs[2])")
    haskey(derivs, m.rhs) ||
        error("maxwell_residual: missing derivative ∂$(m.rhs[1])/∂$(m.rhs[2])")
    return derivs[m.lhs] - m.coeff * derivs[m.rhs]
end
function maxwell_residual(p::ThermodynamicPotential; derivs)
    return maxwell_residual(maxwell_relation(p); derivs=derivs)
end

export MaxwellRelation, maxwell_relation, maxwell_residual
