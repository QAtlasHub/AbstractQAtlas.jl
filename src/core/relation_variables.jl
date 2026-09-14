# core/relation_variables.jl — the RelationVariable layer.
#
# A relation's variables are keyed by the TYPE of the physical thing they
# denote, not by a formula-letter `Symbol` (which drifts and collides — see
# docs/design/type-keyed-interface.md).  The identity-bearing kinds are:
#
#   AbstractQuantity   observables      (core/quantities.jl)
#   AbstractField      control fields   (core/fields.jl)
#   AbstractCoordinate evaluation point (ω, q)              — here
#   AbstractExponent   critical exps    (α, β, γ, ν, η, …)  — here
#
# and a variable also carries a `support` — WHERE it is evaluated — so the same
# type at different regions/points is a distinct key.  The type-keyed prototype
# uses only the trivial `Global` support; `Region` / point / pair supports land
# with the entanglement layer, but the key is `(type, support)` from day one so
# adding them never forces a re-key (design note R1).

"""
    AbstractCoordinate

Parent for evaluation-coordinate variables — a frequency `ω`, a momentum `q`:
the *point at which* a quantity is evaluated, not a subject of the identity.
Coordinates usually appear as lightweight supplied slots rather than typed keys
(design note R3); the type exists so they *can* be keyed when it matters.
"""
abstract type AbstractCoordinate end
export AbstractCoordinate

"""
    Frequency <: AbstractCoordinate

The frequency `ω` at which a dynamical quantity is evaluated.
"""
struct Frequency <: AbstractCoordinate end
export Frequency

"""
    Momentum <: AbstractCoordinate

The momentum / wavevector `q` at which a quantity is evaluated.
"""
struct Momentum <: AbstractCoordinate end
export Momentum

"""
    AbstractExponent

Parent for critical-exponent variables (`α`, `β`, `γ`, `δ`, `ν`, `η`, `z`).
Typing exponents separates the critical-exponent `β` from the inverse
temperature [`InverseTemperature`](@ref) — the two `:β`s a symbol key conflates.
Concrete exponents are introduced when the criticality domain migrates.
"""
abstract type AbstractExponent end
export AbstractExponent
# The concrete exponents the parent above was waiting for. A `VariableKey` is a
# they were. A `VariableKey` is a type, so `β` the order-parameter exponent and
# `β` the inverse temperature were one node in the symbol graph while being two
# quantities: MEASURED, nineteen relations produce `:β` and sixteen of them mean
# the temperature. `α` was likewise shared with the Rényi index and `γ` with the
# topological entanglement entropy. Typing the scaling side separates them without
# touching the others, since only the thirteen `:scaling` relations are annotated.

"""
    SpecificHeatExponent() <: AbstractQuantity

`α`, the specific heat's divergence at a critical point, `c ∼ |t|^{-α}`.
"""
struct SpecificHeatExponent <: AbstractExponent end
export SpecificHeatExponent

"""
    OrderParameterExponent() <: AbstractQuantity

`β`, the order parameter's vanishing, `m ∼ (-t)^β`.  Not
[`InverseTemperature`](@ref), which wears the same letter across most of this
registry.
"""
struct OrderParameterExponent <: AbstractExponent end
export OrderParameterExponent

"""
    SusceptibilityExponent() <: AbstractQuantity

`γ`, the susceptibility's divergence, `χ ∼ |t|^{-γ}`.
"""
struct SusceptibilityExponent <: AbstractExponent end
export SusceptibilityExponent

"""
    CriticalIsothermExponent() <: AbstractQuantity

`δ`, the critical isotherm's shape, `m ∼ h^{1/δ}` at `t = 0`.
"""
struct CriticalIsothermExponent <: AbstractExponent end
export CriticalIsothermExponent

"""
    CorrelationLengthExponent() <: AbstractQuantity

`ν`, the correlation length's divergence, `ξ ∼ |t|^{-ν}`.  The exponent, where
[`CorrelationLength`](@ref) is the length itself.
"""
struct CorrelationLengthExponent <: AbstractExponent end
export CorrelationLengthExponent

"""
    AnomalousDimension() <: AbstractQuantity

`η`, the correlation function's decay at criticality, `G(r) ∼ r^{-(d-2+η)}`.
"""
struct AnomalousDimension <: AbstractExponent end
export AnomalousDimension

"""
    LargeSpinExponent() <: AbstractQuantity

`ζ` of the large-spin fixed point, where the effective moment GROWS under
renormalization ([IgloiMonthus2005](@cite), §A.5).  A random-walk argument on the
signs of the couplings gives `ζ = 1/2`.

Not the correlation-matrix eigenvalue that wears the same letter in
[`EntanglementSpectrumCorrelation`](@ref); one is an exponent and the other an
occupation in `(0, 1)`.
"""
struct LargeSpinExponent <: AbstractExponent end
export LargeSpinExponent

"""
    RelationVariable

Union of the identity-bearing kinds a relation variable can key on:
[`AbstractQuantity`](@ref), [`AbstractField`](@ref), [`AbstractCoordinate`](@ref),
[`AbstractExponent`](@ref).  The *type* is the variable's identity; the formula
letter in the relation body is a private local binding.
"""
const RelationVariable = Union{
    AbstractQuantity,AbstractField,AbstractCoordinate,AbstractExponent
}
export RelationVariable

# Family erasure: `Susceptibility{I}` → `Susceptibility`, so an index-parametric
# quantity collapses to one node/key.  Defined here (the variable-identity file)
# and shared by the auto-derived `quantities` (relations/interface.jl) and the
# quantity graph (structure/graph.jl).
_family(::Type{T}) where {T} = Base.typename(T).wrapper
_family(q::AbstractQuantity) = _family(typeof(q))

# ─── Group slots: how many members of an ABSTRACT group a slot takes ─────
#
# A parametric FAMILY (`Susceptibility`) needs no quantifier: its components are
# one quantity at different indices, so a law written on the family is
# component-agnostic by construction and "check every component" is the only
# reading (§8a).
#
# An abstract GROUP (`AbstractGap`, or a downstream `AbstractTightBindingQuantity`)
# is different in kind — its members are DIFFERENT quantities.  Both readings are
# wanted and they are not interchangeable:
#
#   EachOf{AbstractGap}   every gap is non-negative        — one report row per member
#   AnyOf{AbstractGap}    THE gap sets ξ = 1/Δ             — one member, named by the caller
#
# MEASURED: applying the family rule to a group turns the second kind into false
# violations — `ξ = 1/Δ` on a bag of three gaps reports the two irrelevant ones as
# VIOLATED purely for sharing a supertype.  So the quantifier is written at the
# declaration, and a BARE abstract slot is rejected with a message naming both.

"""
    EachOf{G}

Slot quantifier: the relation holds for **every** member of the abstract group
`G`.  Written as a slot key, `g::EachOf{AbstractGap}`, it behaves exactly like a
parametric-family slot — [`relation_report`](@ref) auto-discovers each concrete
member present in a bag and emits one row per member.

Use it only when the law really is member-agnostic (`every gap ≥ 0`).  For a law
about one member, use [`AnyOf`](@ref) — stating it with `EachOf` reports the
other members as violated purely because they share a supertype.
"""
struct EachOf{G} end
export EachOf

"""
    AnyOf{G}

Slot quantifier: the relation holds for **one** member of the abstract group `G`,
and the caller says which (`check(rel, b; subject = MassGap)`).

Unlike [`EachOf`](@ref), such a relation is *not auto-discoverable*: the engine
cannot know which member a one-member law is about, so [`relation_report`](@ref)
does not instantiate it and [`applicable_relations`](@ref) does not list it.
[`ambiguous_relations`](@ref) lists them instead — the pending work is visible
rather than silently absent.
"""
struct AnyOf{G} end
export AnyOf

# The abstract group a quantifier wraps, or `nothing` for anything else.
_group(::Type{EachOf{G}}) where {G} = G
_group(::Type{AnyOf{G}}) where {G} = G
_group(@nospecialize(T)) = nothing

# ─── Support: WHERE a variable is evaluated ─────────────────────────────

"""
    Support

Where / on what a variable is evaluated.  [`Global`](@ref) — the whole system,
no decoration — is the default and the only support the type-keyed prototype
uses; region / point / pair supports arrive with the entanglement layer.

!!! note "Equality contract"
    A [`VariableKey`](@ref) is a `Dict` key, so it hashes and compares by its
    `support`.  Every `Support` subtype MUST therefore implement value-based
    `Base.==` and `Base.hash` (Julia's struct default is identity `===`).
    [`Global`](@ref) satisfies this trivially as a zero-field singleton; a future
    `Region = Set{AbstractSite}` must define them explicitly, or two
    content-identical regions built separately would key distinct bag entries.
"""
abstract type Support end
export Support

"""
    Global <: Support

The trivial support: a bulk, whole-system quantity with no region/point
decoration.  The default support of every variable.
"""
struct Global <: Support end
export Global
Base.show(io::IO, ::Global) = print(io, "global")

"""
    OrderSupport(order) <: Support

The support of a quantity that is a one-parameter FAMILY — the Rényi entropy `S_α`,
the Tsallis entropy `S_q` — where the order is what distinguishes one member from
another.

It exists because a [`VariableKey`](@ref) is `(type, support)` and the type alone
cannot tell two orders apart. Carrying the order in a plain field does NOT help:
`_as_key` builds the key from `typeof(v)`, and `typeof` erases a non-parametric
struct's field, so before this existed

    bag(RenyiEntropy(2) => 0.5, RenyiEntropy(3) => 0.7)

silently kept only `0.7` — same key, second write wins, no error. MEASURED.

This is the order twin of [`RegionSupport`](@ref): the support slot already existed
for exactly this purpose, "same quantity, different instance", so a one-parameter
family belongs in it rather than in a new type parameter per order.
"""
struct OrderSupport{T} <: Support
    order::T
end
Base.:(==)(a::OrderSupport, b::OrderSupport) = a.order == b.order
Base.hash(a::OrderSupport, h::UInt) = hash(a.order, hash(:OrderSupport, h))
Base.show(io::IO, s::OrderSupport) = print(io, "order ", s.order)
export OrderSupport

"""
    SizeSupport(L) <: Support

The support of a quantity measured on a FINITE SYSTEM of linear size `L`, where
the size is what distinguishes one measurement from another.

The size twin of [`OrderSupport`](@ref), and required for the same reason: a
[`VariableKey`](@ref) is `(type, support)`, so without it a finite-size sweep
cannot even be written down. MEASURED, on the present bag:

    bag(Typical(MassGap()) => 1e-2, Typical(MassGap()) => 1e-4)
    # ERROR: bag: duplicate key VariableKey(Typical{MassGap}). Two entries claim
    # the same identity slot; if they are different quantities, one of them needs
    # a support that says so.

Build the key with [`at_size`](@ref).
"""
struct SizeSupport{T} <: Support
    size::T
end
Base.:(==)(a::SizeSupport, b::SizeSupport) = a.size == b.size
Base.hash(a::SizeSupport, h::UInt) = hash(a.size, hash(:SizeSupport, h))
Base.show(io::IO, s::SizeSupport) = print(io, "L = ", s.size)
export SizeSupport

"""
    at_size(q::AbstractQuantity, L) -> VariableKey
    at_size(::Type{<:AbstractQuantity}, L) -> VariableKey

The bag key for `q` measured on a system of size `L`, so one bag can hold a whole
finite-size sweep:

```julia
b = bag(at_size(Typical(MassGap()), 16) => 1e-2,
        at_size(Typical(MassGap()), 32) => 1e-4)
```

The twin of [`entanglement_entropy`](@ref) for the size axis: it writes the key
directly, because the size belongs to the measurement and not to the quantity.
A quantity that already needs a support of its own is refused rather than
silently losing it.
"""
function at_size(q::AbstractQuantity, L)
    L isa Real && isfinite(L) ||
        error("at_size: a size must be a finite real, got $L::$(typeof(L)).")
    # `typeof` erases a field, so a quantity whose identity lives in one would key
    # two different measurements to one slot. `RenyiEntropy(2)` and `RenyiEntropy(3)`
    # at the same size collide loudly, and at different sizes fuse into a sweep that
    # never existed. Supports do not compose here, so the pair is refused.
    variable_support(q) isa Global || error(
        "at_size: $(typeof(q)) already keys under $(variable_support(q)), and a key " *
        "carries one support. Size and $(nameof(typeof(variable_support(q)))) cannot " *
        "be combined, so this measurement has no slot to go in.",
    )
    return VariableKey(typeof(q), SizeSupport(L))
end
function at_size(@nospecialize(Q::Type), L)
    L isa Real && isfinite(L) ||
        error("at_size: a size must be a finite real, got $L::$(typeof(L)).")
    # The instance method's guard applies here too: a type whose instances key under
    # a support of their own would lose it just as silently through this spelling.
    Q <: AbstractQuantity &&
        !(variable_support(Q()) isa Global) &&
        error(
            "at_size: $Q keys under $(variable_support(Q())), and a key carries one " *
            "support, so this measurement has no slot to go in.",
        )
    return VariableKey(Q, SizeSupport(L))
end
export at_size

"""
    variable_support(v) -> Support

The support a variable INSTANCE keys under. `Global()` for everything whose type is
its whole identity; overridden where an instance carries data that distinguishes it
from another instance of the same type.

Add an override here whenever a quantity gains such a field — that is the one place
that decides whether two instances share a bag slot, and the failure mode of getting
it wrong is silent (see [`OrderSupport`](@ref)).
"""
variable_support(::RelationVariable) = Global()
variable_support(q::RenyiEntropy) = OrderSupport(q.α)
variable_support(q::TsallisEntropy) = OrderSupport(q.q)
export variable_support

# ─── VariableKey: (type, support) — the collision-proof identity ────────

"""
    VariableKey(type::Type, support::Support = Global())

The identity of a relation variable: the quantity / field / coordinate /
exponent `type` together with its `support`.  This — never a formula-letter
`Symbol` — is what the type-keyed bag, [`relation_report`](@ref), and the
derivation graph match on, so distinct types (and, later, distinct supports of
one type) can never collide.
"""
struct VariableKey
    type::Type{<:RelationVariable}   # constrained: a typo'd non-variable type errors here
    support::Support
end
VariableKey(t::Type) = VariableKey(t, Global())

Base.:(==)(a::VariableKey, b::VariableKey) = a.type === b.type && a.support == b.support
Base.hash(k::VariableKey, h::UInt) = hash(k.type, hash(k.support, hash(:VariableKey, h)))
function Base.show(io::IO, k::VariableKey)
    # print the FULL type (with parameters) — `Susceptibility{(:z, :z)}`, never a
    # parameter-dropping `Susceptibility`, so two components are never confused.
    print(io, "VariableKey(", k.type)
    k.support isa Global || print(io, ", ", k.support)
    return print(io, ")")
end
export VariableKey
