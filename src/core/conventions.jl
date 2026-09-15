# core/conventions.jl: the convention layer.
#
# A `VariableKey` pins WHICH quantity and WHERE it is evaluated. It does not pin
# how the number is written, so one type at one support is still two different
# values when one calculation counts bits and the other nats. Today that is
# prose: relations/entanglement.jl's header says "Entropies here are in NATS"
# and nothing reads it.
#
# Three parts, kept separate because they are declared by different people:
# the owner of a quantity says what the relations are written in
# (`canonical_convention`), a project says what ITS numbers are written in
# (`conventions`), and `bag` converts between them on entry so an undeclared
# value can never reach a relation.
#
# Conversion is opt-in per convention pair and refuses otherwise, because most
# pairs are not a rescaling: Pauli and spin-1/2 operators move each term of a
# Hamiltonian by a different factor, and a disorder average is not a rescaled
# typical value.

"""
    Convention

Parent for the ways one quantity's value can be written: the base of the
logarithm in an entropy, the normalisation of the operators an energy is built
from, the statistic a disorder average reports.

Open by design. A new axis is a new subtype plus a [`convert_convention`](@ref)
method and needs no change to this file, which is why the package ships the
protocol rather than a list of axes. A quantity with two independent axes gets
one subtype carrying both, since [`canonical_convention`](@ref) answers with a
single value.
"""
abstract type Convention end
export Convention

"""
    canonical_convention(Q::Type) -> Union{Convention,Nothing}

The convention this package's relations are written in for `Q`, or `nothing`
when `Q` has no convention axis at all.

`nothing` is the default and is deliberately not filled in by supertype: the
Tsallis entropy is `(1 - Tr ρ^q)/(q - 1)`, which carries no logarithm, so
declaring a base for every [`AbstractEntanglementMeasure`](@ref) would give it
an axis it does not have. Same opt-in reasoning as
[`obeys_entropy_inequalities`](@ref).
"""
canonical_convention(::Type{<:RelationVariable}) = nothing
export canonical_convention

"""
    convert_convention(to::Convention, from::Convention, Q::Type, v)

`v`, written in `from`, expressed in `to`.

Equal conventions return `v` untouched, so an exact-arithmetic value stays
exact. Anything else refuses unless a method says the pair converts: two
conventions that are not a rescaling have no conversion, and passing `v`
through there would hand a relation a number from the wrong statistic.
"""
function convert_convention(to::Convention, from::Convention, @nospecialize(Q::Type), v)
    to == from && return v
    return error(
        "convert_convention: no conversion from $from to $to for $Q. Define a " *
        "`convert_convention` method if the two are a rescaling of each other; " *
        "if they are not, the values have to be recomputed rather than converted.",
    )
end
export convert_convention

# ─── The log-base axis ───────────────────────────────────────────────────

"""
    LogBase(base::Real) <: Convention

The base of the logarithm an entropy is measured with: [`Nats`](@ref) is
`LogBase(ℯ)` and [`Bits`](@ref) is `LogBase(2)`.

A coefficient multiplying a logarithm (`c`, `c̃`) is unchanged by the base, since
it rescales entropy and logarithm alike. An additive constant (`c₁`, `ln g`) and
a bare difference of entropies are not, and that is where a transcribed formula
silently gains or loses a factor of `ln 2`.
"""
struct LogBase <: Convention
    base::Float64
    function LogBase(base::Real)
        base > 0 && base != 1 ||
            throw(ArgumentError("LogBase: base must be positive and not 1; got $base"))
        return new(Float64(base))
    end
end
export LogBase

"""
    Nats

`LogBase(ℯ)`, the convention every entropy relation in this package is written in.
"""
const Nats = LogBase(ℯ)
export Nats

"""
    Bits

`LogBase(2)`, `S = -Tr ρ log₂ ρ`, which is what the entanglement literature
usually counts.
"""
const Bits = LogBase(2)
export Bits

function convert_convention(to::LogBase, from::LogBase, @nospecialize(Q::Type), v)
    to == from && return v
    return v * (log(from.base) / log(to.base))
end

# The twelve entanglement measures whose ABQ definition contains a logarithm, or
# is an additive combination of ones that do. The four that are absent are absent
# because their defining relation has no logarithm to take a base of:
# `TsallisEntropy` is `(1 - Tr ρ^q)/(q - 1)`, `Concurrence` is a wavefunction
# amplitude, and `Tangle`/`ThreeTangle` are built from it by squaring.
for Q in (
    :VonNeumannEntropy,
    :FermionicEntanglementEntropy,
    :RenyiEntropy,
    :MutualInformation,
    :ConditionalEntropy,
    :RelativeEntropy,
    :MeasurementEntropy,
    :MarkovEntropy,
    :TripartiteInformation,
    :TopologicalEntanglementEntropy,
    :LogarithmicNegativity,
    :PageEntropy,
)
    @eval canonical_convention(::Type{$Q}) = Nats
end

# ─── What a project declares ─────────────────────────────────────────────

"""
    ConventionSet

What one project's numbers are written in, as a map from quantity type to
[`Convention`](@ref). Build one with [`conventions`](@ref) and hand it to
[`bag`](@ref).

Declared once next to the calculation, not repeated at each call, because the
convention is a property of how the numbers were produced.
"""
struct ConventionSet
    declared::Dict{Type,Convention}
end
export ConventionSet

"""
    conventions(pairs...) -> ConventionSet

Declare what this project's values are written in:

```julia
conventions(VonNeumannEntropy => Bits)              # this one quantity
conventions(AbstractEntanglementMeasure => Bits)    # every entropy that has a base
```

A key may be a concrete type, an INSTANCE (reduced to its type), or an abstract
supertype. The supertype form is the usable one when a bag holds several
entropies: it reaches every subtype that declares a
[`canonical_convention`](@ref) and skips the ones that have no such axis, so
naming `AbstractEntanglementMeasure` does not claim a base for the Tsallis
entropy. Naming a concrete type that has no axis is an error, since that is a
claim about that type.

Lookup is most-specific-first, so a concrete entry overrides a supertype entry.
"""
function conventions(pairs::Pair...)
    d = Dict{Type,Convention}()
    for (k, c) in pairs
        T = k isa Type ? k : typeof(k)
        T <: RelationVariable ||
            error("conventions: $T is not a relation variable, so no relation reads it.")
        c isa Convention ||
            error("conventions: the value for $T is $(typeof(c)), not a Convention.")
        haskey(d, T) && error(
            "conventions: duplicate key $T. A quantity is written in one " *
            "convention; two entries for it cannot both be what the numbers are.",
        )
        # A concrete type with no axis is a claim about that type, so it is
        # refused here rather than silently ignored at conversion time.
        if isconcretetype(T) && canonical_convention(T) === nothing
            error(
                "conventions: $T declares no convention axis, so $c cannot be " *
                "converted away from. If $T does have one, give it a " *
                "`canonical_convention` method; if the declaration was meant for " *
                "its siblings, key it on the supertype instead.",
            )
        end
        d[T] = c
    end
    return ConventionSet(d)
end
export conventions

"""
    declared_convention(cs::ConventionSet, Q::Type) -> Union{Convention,Nothing}

What `cs` says `Q`'s values are written in, taking the most specific entry that
`Q` is a subtype of; `nothing` when nothing in `cs` covers `Q`.

Matched by `<:`, not by walking `supertype`, because a parametric quantity's
supertype chain SKIPS its own family: `supertype(Energy{:per_site})` is
`AbstractThermalPotential`, so a walk never reaches the `Energy` a project keyed
its declaration on, and the value goes into the bag unconverted.

Two declared types that both cover `Q` and are unrelated to each other are
refused rather than resolved by `Dict` order.
"""
function declared_convention(cs::ConventionSet, @nospecialize(Q::Type))
    best, bestT = nothing, nothing
    for (T, c) in cs.declared
        Q <: T || continue
        if bestT === nothing || T <: bestT
            best, bestT = c, T
        elseif !(bestT <: T)
            error(
                "declared_convention: $Q is covered by both $bestT and $T, which are " *
                "unrelated, so neither is the more specific. Key the declaration on " *
                "whichever one the values were actually written in.",
            )
        end
    end
    return best
end
export declared_convention

"""
    in_canonical_convention(cs::ConventionSet, Q::Type, v)

`v` rewritten in the convention this package's relations use for `Q`.

Returns `v` untouched when `cs` says nothing about `Q`, or when `Q` has no
convention axis, which is what makes the layer opt-in: a project that declares
nothing gets the behaviour it had before.
"""
function in_canonical_convention(cs::ConventionSet, @nospecialize(Q::Type), v)
    from = declared_convention(cs, Q)
    from === nothing && return v
    to = canonical_convention(Q)
    to === nothing && return v
    return convert_convention(to, from, Q, v)
end
export in_canonical_convention
