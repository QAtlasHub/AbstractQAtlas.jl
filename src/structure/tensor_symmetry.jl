# structure/tensor_symmetry.jl — symmetries of the response tensors.
#
# A nonlinear susceptibility is *essentially* a higher-order tensor, and
# its higher-order structure carries a known, model-independent symmetry:
# **intrinsic permutation symmetry**.  The order-n response
# `χ⁽ⁿ⁾_{α; β₁…βₙ}(ω₁, …, ωₙ)` is invariant under any simultaneous
# permutation of the (field-index, frequency) pairs
# `{(β₁, ω₁), …, (βₙ, ωₙ)}` — a consequence of the fields being
# interchangeable perturbations (Armstrong, Bloembergen, Ducuing &
# Pershan, Phys. Rev. 127, 1918 (1962)).  This layer makes that symmetry
# explicit at the type level, so a consumer knows which components are
# forced equal without re-deriving it.  Because the symmetry pairs each
# field index with its frequency, `field_permutation` reports the
# permutation the symmetry also applies to the frequency arguments — the
# piece a frequency-resolved (nonlinear/dynamical) verification actually
# needs, and which an index-only canonicalization silently drops.

"""
    intrinsic_permutation_symmetric(quantity) -> Bool
    intrinsic_permutation_symmetric(::Type{<:AbstractQuantity}) -> Bool

Whether the response tensor is invariant under permutation of its field
indices (paired with their frequencies) — the intrinsic permutation
symmetry of nonlinear response.  `true` for the susceptibilities and the
conductivity, `false` otherwise (the default).
"""
function intrinsic_permutation_symmetric(q::AbstractQuantity)
    return intrinsic_permutation_symmetric(typeof(q))
end
intrinsic_permutation_symmetric(::Type{<:AbstractQuantity}) = false
intrinsic_permutation_symmetric(::Type{<:Susceptibility}) = true
intrinsic_permutation_symmetric(::Type{<:DynamicalSusceptibility}) = true
intrinsic_permutation_symmetric(::Type{<:Conductivity}) = true
intrinsic_permutation_symmetric(::Type{<:DynamicalConductivity}) = true
export intrinsic_permutation_symmetric

"""
    canonical_component(χ) -> typeof(χ)

The canonical representative of `χ`'s intrinsic-permutation-symmetry
class: the same response tensor with its **field** indices sorted (the
response index — the first — is fixed).  Two components are forced equal
by the symmetry iff they share a `canonical_component`
([`permutation_equivalent`](@ref)):

```julia
canonical_component(Susceptibility(:x, :z, :y)) === Susceptibility(:x, :y, :z)
```
"""
function canonical_component(χ::Susceptibility)
    i = indices(χ)
    return Susceptibility(i[1], sort!(collect(i[2:end]))...)
end
function canonical_component(χ::DynamicalSusceptibility)
    i = indices(χ)
    return DynamicalSusceptibility(i[1], sort!(collect(i[2:end]))...)
end
function canonical_component(χ::Conductivity)
    i = indices(χ)
    return Conductivity(i[1], sort!(collect(i[2:end]))...)
end
function canonical_component(χ::DynamicalConductivity)
    i = indices(χ)
    return DynamicalConductivity(i[1], sort!(collect(i[2:end]))...)
end
export canonical_component

"""
    permutation_equivalent(a, b) -> Bool

Whether response tensors `a` and `b` are forced equal by intrinsic
permutation symmetry — i.e. differ only by a permutation of their field
indices.  `χ⁽²⁾_{x;yz} = χ⁽²⁾_{x;zy}`, so they are equivalent; `χ_{x;yz}`
and `χ_{y;xz}` are not (different response index).

!!! warning "the symmetry pairs field indices with frequencies"
    The intrinsic permutation symmetry acts on **(field-index, frequency)
    pairs**, so for a frequency-resolved response the equality holds only
    when the frequency arguments are permuted to match:
    `χ⁽²⁾_{x;yz}(ω₁, ω₂) == χ⁽²⁾_{x;zy}(ω₂, ω₁)`, **not**
    `χ⁽²⁾_{x;zy}(ω₁, ω₂)`.  A static [`Susceptibility`](@ref)
    (`frequency_arguments == 0`) has no frequencies to permute, so
    permutation-equivalent components are numerically equal outright; for
    a [`DynamicalSusceptibility`](@ref) / [`Conductivity`](@ref)
    (`frequency_arguments > 0`) apply [`field_permutation`](@ref) to the
    frequency arguments before comparing.
"""
function permutation_equivalent(a::AbstractQuantity, b::AbstractQuantity)
    (intrinsic_permutation_symmetric(a) && typeof(a).name === typeof(b).name) ||
        return false
    return canonical_component(a) === canonical_component(b)
end
export permutation_equivalent

"""
    field_permutation(χ) -> NTuple{n,Int}

The permutation `π` of `χ`'s `n = response_order(χ)` field slots that
brings its field indices to `canonical_component` (sorted) order —
`Tuple(sortperm(collect(indices(χ)[2:end])))`.

Because the intrinsic permutation symmetry pairs each field index with
its frequency, `π` is *also* the permutation the symmetry applies to the
frequency arguments:

    χ_{α; β₁…βₙ}(ω₁, …, ωₙ) == canonical_component(χ)(ω_{π₁}, …, ω_{πₙ})

so it is exactly what a consumer needs to check one frequency-resolved
component against another (or against the canonical representative).  For
a static [`Susceptibility`](@ref) (`frequency_arguments == 0`) there is
nothing to permute and permutation-equivalent components are equal
outright; `π` still reports the field-index sort.

```julia
field_permutation(DynamicalSusceptibility(:x, :z, :y))  # (2, 1) — swap the two frequencies
field_permutation(DynamicalSusceptibility(:x, :y, :z))  # (1, 2) — already canonical
```
"""
function field_permutation(χ::AbstractQuantity)
    intrinsic_permutation_symmetric(χ) || return error(
        "field_permutation: $(typeof(χ)) has no intrinsic permutation symmetry"
    )
    return Tuple(sortperm(collect(indices(χ)[2:end])))
end
export field_permutation
