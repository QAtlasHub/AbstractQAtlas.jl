# structure/spectral.jl — the dynamical-quantity inter-relationship graph.
#
# The frequency-resolved quantities form a chain of *operations*:
#
#   SelfEnergy  --dyson-->  RetardedGreensFunction
#               --neg_im_over_pi-->  SpectralFunction
#               --bz_average-->  DensityOfStates
#
#   DynamicalCorrelation{I}  --spacetime_fourier-->  DynamicalStructureFactor  (linear, 2-point)
#                            --kubo-->  DynamicalSusceptibility{I}             (order-preserving)
#                                       --low_frequency_limit-->  NMRSpinRelaxationRate
#
# The DynamicalCorrelation is the root of the response branch: both the
# symmetric structure factor S(q,ω) and the antisymmetric (commutator)
# response χ(q,ω) come from it — the two sides of the fluctuation–
# dissipation theorem.  It is order-parametric: an n-th order response
# χ⁽ⁿ⁾(ω₁…ωₙ) is the retarded part of the n-time (n+1-point)
# DynamicalCorrelation of the SAME order, so the :kubo edge preserves the
# response order (the structure factor is the linear, 2-point branch).
#
# Each edge records WHICH operation turns the source quantity into the
# target — the inter-quantity genealogy the response tree
# (structure/response.jl) is for thermodynamics.  Where the operation is
# a single (q, ω)-point identity (`neg_im_over_pi`, `dyson`) the edge
# points at the exact `@relation`; where it is a transform / BZ sum /
# limit (`spacetime_fourier`, `bz_average`, `low_frequency_limit`) there
# is no single-point relation — evaluating it is the functional sibling's
# job (issue #14) — so `origin_relation` returns `nothing`.

"""
    SpectralOrigin(from, via)

One edge of the dynamical-quantity graph: the quantity carrying it is
obtained from quantity type `from` by the operation named `via` — one of
`:dyson`, `:neg_im_over_pi`, `:bz_average`, `:spacetime_fourier`,
`:low_frequency_limit`, `:kubo`.
"""
struct SpectralOrigin
    from::Type
    via::Symbol
end

"""
    spectral_origin(quantity) -> Union{SpectralOrigin,Nothing}
    spectral_origin(::Type{<:AbstractQuantity}) -> Union{SpectralOrigin,Nothing}

The dynamical-graph edge of `quantity`: the (`from`, `via`) it is
obtained from, or `nothing` for a source quantity
([`SelfEnergy`](@ref), [`DynamicalCorrelation`](@ref)) or a quantity
outside the graph.

```julia
spectral_origin(DensityOfStates())      # SpectralOrigin(SpectralFunction, :bz_average)
spectral_origin(SpectralFunction())     # SpectralOrigin(RetardedGreensFunction, :neg_im_over_pi)
spectral_origin(RetardedGreensFunction()) # SpectralOrigin(SelfEnergy, :dyson)
```
"""
spectral_origin(q::AbstractQuantity) = spectral_origin(typeof(q))
spectral_origin(::Type{<:AbstractQuantity}) = nothing

spectral_origin(::Type{RetardedGreensFunction}) = SpectralOrigin(SelfEnergy, :dyson)
function spectral_origin(::Type{SpectralFunction})
    return SpectralOrigin(RetardedGreensFunction, :neg_im_over_pi)
end
spectral_origin(::Type{DensityOfStates}) = SpectralOrigin(SpectralFunction, :bz_average)
function spectral_origin(::Type{DynamicalStructureFactor})
    return SpectralOrigin(DynamicalCorrelation, :spacetime_fourier)
end
# Kubo formula: the dynamical susceptibility is the retarded (commutator)
# response function of the correlation.  Linear χ(ω) = FT of
# (i/ħ)θ(t)⟨[A(t),B(0)]⟩ (Kubo, J. Phys. Soc. Jpn. 12, 570 (1957)); the
# n-th order χ⁽ⁿ⁾(ω₁…ωₙ) is the multi-time n-fold nested-commutator
# response (2D coherent spectroscopy — Wan & Armitage, PRL 122, 257401
# (2019)).  An n-th order response is an n-TIME correlation: order-n
# χ⁽ⁿ⁾{I} routes to the same-order DynamicalCorrelation{I} (the (n+1)-point
# kernel), so `frequency_arguments` agree on both sides of the edge.
function spectral_origin(::Type{DynamicalSusceptibility{I}}) where {I}
    return SpectralOrigin(DynamicalCorrelation{I}, :kubo)
end
# the index-erased family (e.g. reached via NMRSpinRelaxationRate) routes to
# the correlation family
function spectral_origin(::Type{DynamicalSusceptibility})
    return SpectralOrigin(DynamicalCorrelation, :kubo)
end
# Kubo edge for the AC conductivity: σ⁽ⁿ⁾{I} is the retarded part of the
# same-order n-time current–current correlation (order-faithful, like the
# susceptibility) — Kubo, J. Phys. Soc. Jpn. 12, 570 (1957).
function spectral_origin(::Type{DynamicalConductivity{I}}) where {I}
    return SpectralOrigin(CurrentCorrelation{I}, :kubo)
end
spectral_origin(::Type{DynamicalConductivity}) = SpectralOrigin(CurrentCorrelation, :kubo)
# the current-noise spectral density is the space-time FT of the current
# correlation — the fluctuation side, mirroring DynamicalStructureFactor
function spectral_origin(::Type{<:CurrentNoise})
    return SpectralOrigin(CurrentCorrelation, :spacetime_fourier)
end
function spectral_origin(::Type{NMRSpinRelaxationRate})
    return SpectralOrigin(DynamicalSusceptibility, :low_frequency_limit)
end
export spectral_origin, SpectralOrigin

"""
    origin_relation(via::Symbol) -> Union{AbstractRelation,Nothing}

The exact single-`(q, ω)`-point [`@relation`](@ref) that realizes the
operation `via`, or `nothing` when the operation is a transform / sum /
limit with no pointwise form (its evaluation belongs to the functional
sibling, issue #14):

- `:dyson` → [`Dyson`](@ref),
- `:neg_im_over_pi` → [`SpectralFromGreens`](@ref),
- `:bz_average`, `:spacetime_fourier`, `:low_frequency_limit`, `:kubo`
  → `nothing` (transform / sum / limit / commutator-response — no
  single-point form; evaluation is the functional sibling's job).

The Kubo edge (`:kubo`) is a transform of a multi-time correlation
(Kubo, J. Phys. Soc. Jpn. 12, 570 (1957)), so it has no single-`(q,ω)`-
point relation here.
"""
origin_relation(via::Symbol) = origin_relation(Val(via))
origin_relation(::Val{:dyson}) = Dyson()
origin_relation(::Val{:neg_im_over_pi}) = SpectralFromGreens()
origin_relation(::Val) = nothing
export origin_relation

"""
    operation_scope(via::Symbol) -> Symbol

Which layer owns the dynamical-graph operation `via` — the **scope line**
of issue #14:

- `:definitional` — a **pointwise** identity relating quantity *values* at
  a single `(q, ω)` (or a supplied scalar: an integral, a derivative). It
  lives HERE, in this stdlib-only definitional package, as an
  [`@relation`](@ref) (`:dyson`, `:neg_im_over_pi`; and every
  supplied-integral / supplied-derivative relation such as the sum rules
  and Kramers–Kronig).
- `:functional` — a **transform / sum / limit** that must represent a
  quantity as a FUNCTION and act on it globally (a BZ average, a
  space-time Fourier transform, an ω → 0 limit, the Kubo response). Its
  *evaluation* belongs to the future ParaLA-based functional sibling;
  only its structural edge lives here (`spectral_origin`).

The line is exactly `origin_relation`'s split: an operation is
`:definitional` iff it has a pointwise `@relation`.  Grey zone (issue #14,
cf. #6): a sum rule is `:definitional` — the RELATION checks a supplied
number here, while COMPUTING that number from the function is
`:functional` (the sibling's job).
"""
function operation_scope(via::Symbol)
    return origin_relation(via) === nothing ? :functional : :definitional
end
export operation_scope

"""
    spectral_chain(quantity) -> Vector{Any}

The dynamical-graph path from `quantity` back to its source, as the list
of quantity types `[typeof(quantity), from, …, source]`.  A source (or
off-graph) quantity returns the singleton `[typeof(quantity)]`.

```julia
spectral_chain(DensityOfStates())
# [DensityOfStates, SpectralFunction, RetardedGreensFunction, SelfEnergy]
# i.e. ρ ⟵ A ⟵ G^R ⟵ Σ : the density of states is built from the
# self-energy through Dyson, the spectral representation, and the BZ sum.
```
"""
function spectral_chain(q::AbstractQuantity)
    chain = Any[typeof(q)]
    cur = typeof(q)
    for _ in 1:32   # depth cap guards a mis-declared cycle
        o = spectral_origin(cur)
        o === nothing && return chain
        push!(chain, o.from)
        cur = o.from
    end
    return error("spectral_chain: graph exceeded depth cap — cyclic edge?")
end
export spectral_chain
