# structure/fourier.jl — which representation each quantity lives in, and
# the Fourier-conjugate quantity it maps to.
#
# `representation(q)` tags a quantity with its spatial and/or temporal
# space; `fourier_conjugate_quantity(q)` names the quantity obtained by
# Fourier-transforming ALL of `q`'s representations to their conjugates
# (the momentum/frequency object ↔ the real-space/time object).  The
# actual (discrete) transform is `AbstractFFTs.fft`/`ifft` on a grid —
# the functional sibling's job (issue #14); here it is structure.

"""
    representation(quantity) -> Tuple{Vararg{AbstractRepresentation}}
    representation(::Type{<:AbstractQuantity}) -> Tuple

The spatial and/or temporal representation(s) a quantity is expressed in
— e.g. `(MomentumSpace(), FrequencyDomain())` for `S(q, ω)`,
`(RealSpace(),)` for a real-space correlation, `()` for a global
thermodynamic quantity with no space/time resolution.
"""
representation(q::AbstractQuantity) = representation(typeof(q))
representation(::Type{<:AbstractQuantity}) = ()

# real-space vs momentum-space correlation/structure-factor families
representation(::Type{<:AbstractTwoPointCorrelation}) = (RealSpace(),)
representation(::Type{<:AbstractStructureFactor}) = (MomentumSpace(),)
# the dynamical (space+time / momentum+frequency) family
representation(::Type{<:DynamicalCorrelation}) = (RealSpace(), TimeDomain())
representation(::Type{DynamicalStructureFactor}) = (MomentumSpace(), FrequencyDomain())
function representation(::Type{DynamicalSusceptibility{I}}) where {I}
    return (MomentumSpace(), FrequencyDomain())
end
# the AC conductivity lives in (q, ω); its current–current correlation in (r, t)
representation(::Type{<:DynamicalConductivity}) = (MomentumSpace(), FrequencyDomain())
representation(::Type{<:CurrentCorrelation}) = (RealSpace(), TimeDomain())
# the current-noise spectral density is the (q, ω) FT of the current correlation
representation(::Type{<:CurrentNoise}) = (MomentumSpace(), FrequencyDomain())
# the nonlinear response kernel is the time-domain side of the susceptibility/conductivity
function representation(::Type{ResponseKernel{I}}) where {I}
    return (RealSpace(), TimeDomain())
end
function representation(::Type{CurrentResponseKernel{I}}) where {I}
    return (RealSpace(), TimeDomain())
end
# single-particle propagators live in (q, ω)
representation(::Type{RetardedGreensFunction}) = (MomentumSpace(), FrequencyDomain())
representation(::Type{SelfEnergy}) = (MomentumSpace(), FrequencyDomain())
representation(::Type{SpectralFunction}) = (MomentumSpace(), FrequencyDomain())
export representation

"""
    fourier_conjugate_quantity(quantity) -> Type

The quantity obtained by Fourier-transforming `quantity` in every
representation it carries — the real-space/time object of a
momentum/frequency one and vice versa:

```julia
fourier_conjugate_quantity(StaticStructureFactor)   # SpinCorrelation   (spatial FT)
fourier_conjugate_quantity(DynamicalStructureFactor)# DynamicalCorrelation (space-time FT)
```

Its `representation` is the tuple of `fourier_conjugate`s of the
original's.  Defined for the quantities with an unambiguous conjugate
partner; others (representation-agnostic ones) have none.
"""
fourier_conjugate_quantity(q::AbstractQuantity) = fourier_conjugate_quantity(typeof(q))
fourier_conjugate_quantity(::Type{<:AbstractQuantity}) = nothing

fourier_conjugate_quantity(::Type{StaticStructureFactor}) = SpinCorrelation
fourier_conjugate_quantity(::Type{<:SpinCorrelation}) = StaticStructureFactor
fourier_conjugate_quantity(::Type{DynamicalStructureFactor}) = DynamicalCorrelation
fourier_conjugate_quantity(::Type{<:DynamicalCorrelation}) = DynamicalStructureFactor
# current channel: S^j(q,ω) ↔ ⟨jj⟩(r,t)
fourier_conjugate_quantity(::Type{<:CurrentNoise}) = CurrentCorrelation
fourier_conjugate_quantity(::Type{<:CurrentCorrelation}) = CurrentNoise
# the RESPONSE channel had no Fourier partner at all: the time-domain kernel a real-time
# method computes was unnameable. Each pair is declared FOUR times, not two: both directions
# (Fourier conjugacy is symmetric), and for each direction both the applied type AND the
# index-erased family — exactly as the Kubo edges in structure/spectral.jl are. The family
# form is not decoration. The reflection sweep in test/core/test_invariants.jl reaches leaves
# as bare `UnionAll`s, and a `::Type{X{I}} where {I}` method does not match `f(X)`, so a
# family declaring only the applied form is SKIPPED by the very test meant to catch a
# one-sided declaration.
fourier_conjugate_quantity(::Type{DynamicalSusceptibility{I}}) where {I} = ResponseKernel{I}
fourier_conjugate_quantity(::Type{DynamicalSusceptibility}) = ResponseKernel
fourier_conjugate_quantity(::Type{ResponseKernel{I}}) where {I} = DynamicalSusceptibility{I}
fourier_conjugate_quantity(::Type{ResponseKernel}) = DynamicalSusceptibility
function fourier_conjugate_quantity(::Type{DynamicalConductivity{I}}) where {I}
    return CurrentResponseKernel{I}
end
fourier_conjugate_quantity(::Type{DynamicalConductivity}) = CurrentResponseKernel
function fourier_conjugate_quantity(::Type{CurrentResponseKernel{I}}) where {I}
    return DynamicalConductivity{I}
end
fourier_conjugate_quantity(::Type{CurrentResponseKernel}) = DynamicalConductivity
export fourier_conjugate_quantity

"""
    fourier_pair(a, b) -> Bool

Whether quantities `a` and `b` are Fourier conjugates — the same physics
in conjugate representations (their `representation`s are elementwise
[`fourier_conjugate`](@ref)s, and one is the declared
[`fourier_conjugate_quantity`](@ref) of the other).

```julia
fourier_pair(StaticStructureFactor(), SpinCorrelation(:z, :z))       # true (S(q) ↔ ⟨SS⟩(r))
fourier_pair(DynamicalStructureFactor(), DynamicalCorrelation(:x, :x)) # true (space-time FT)
```
"""
function fourier_pair(a::AbstractQuantity, b::AbstractQuantity)
    ca = fourier_conjugate_quantity(a)
    ca === nothing && return false
    # `b`'s family matches the declared conjugate (allow parametric families),
    # and their representations are elementwise conjugate.
    typeof(b) <: ca || return false
    ra, rb = representation(a), representation(b)
    length(ra) == length(rb) || return false
    return all(fourier_conjugate(x) === y for (x, y) in zip(ra, rb))
end
export fourier_pair
