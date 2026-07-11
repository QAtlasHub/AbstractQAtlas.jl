# structure/response.jl — the response-function derivative genealogy.
#
# Every response function is a derivative of another quantity with
# respect to a control field, and the chain of such derivatives roots at
# the free energy (equivalently, at the partition function via
# `F = −β⁻¹ ln Z`).  This layer encodes that genealogy as a graph of
# edges "child = ∂(parent)/∂(field)".
#
# It records the derivative TOPOLOGY only — which quantity differentiates
# into which, with respect to what.  The exact prefactors and signs live
# in the corresponding `relations/` (the free-energy Legendre / Gibbs–
# Helmholtz / FDT identities): definition here, exact formula there, the
# same split the rest of the package keeps.
#
# The genealogy is *why* the susceptibility is a second field-derivative
# of the free energy: `χ = ∂M/∂h` and `M = −∂F/∂h`, so tracing
# `Susceptibility → Magnetization → FreeEnergy` gives `∂²F/∂h²` — and, in
# the critical region, the same chain relates the exponents β, γ, δ.

"""
    DerivativeEdge(parent, field)

One edge of the response genealogy: the quantity carrying this edge is,
up to a model-independent prefactor, `∂(parent)/∂(field)` — a derivative
of quantity type `parent` with respect to field type `field`.  The exact
prefactor/sign is supplied by the corresponding relation (e.g.
[`GibbsHelmholtz`](@ref) for the energy's `β`-edge).
"""
struct DerivativeEdge
    parent::Type
    field::Type
end

"""
    derivative_edge(quantity) -> Union{DerivativeEdge,Nothing}
    derivative_edge(::Type{<:AbstractQuantity}) -> Union{DerivativeEdge,Nothing}

The genealogy edge of `quantity`: the (`parent`, `field`) it is a
derivative of, or `nothing` for a root potential (the free energy) or a
quantity outside the thermodynamic-derivative tree.

```julia
derivative_edge(SusceptibilityZZ())  # DerivativeEdge(MagnetizationZ, MagneticField)  (χ = ∂M/∂h)
derivative_edge(MagnetizationZ())    # DerivativeEdge(FreeEnergy, MagneticField)      (M = −∂F/∂h)
derivative_edge(FreeEnergy())        # nothing — the root
```
"""
derivative_edge(q::AbstractQuantity) = derivative_edge(typeof(q))
derivative_edge(::Type{<:AbstractQuantity}) = nothing

# Field-derivative branch: M = −∂F/∂h, χ = ∂M/∂h  (per axis component)
derivative_edge(::Type{MagnetizationX}) = DerivativeEdge(FreeEnergy, MagneticField)
derivative_edge(::Type{MagnetizationY}) = DerivativeEdge(FreeEnergy, MagneticField)
derivative_edge(::Type{MagnetizationZ}) = DerivativeEdge(FreeEnergy, MagneticField)
derivative_edge(::Type{SusceptibilityXX}) = DerivativeEdge(MagnetizationX, MagneticField)
derivative_edge(::Type{SusceptibilityYY}) = DerivativeEdge(MagnetizationY, MagneticField)
derivative_edge(::Type{SusceptibilityZZ}) = DerivativeEdge(MagnetizationZ, MagneticField)
# Temperature branch: S = −∂F/∂T, U = ∂(βF)/∂β, C = ∂U/∂T
derivative_edge(::Type{ThermalEntropy}) = DerivativeEdge(FreeEnergy, Temperature)
derivative_edge(::Type{<:Energy}) = DerivativeEdge(FreeEnergy, InverseTemperature)
derivative_edge(::Type{SpecificHeat}) = DerivativeEdge(Energy, Temperature)
export derivative_edge, DerivativeEdge

"""
    is_response(quantity) -> Bool

Whether `quantity` is a response function — i.e. a derivative of another
quantity in the genealogy ([`derivative_edge`](@ref) is non-`nothing`).
"""
is_response(q::AbstractQuantity) = derivative_edge(q) !== nothing
export is_response

"""
    differentiation_chain(quantity) -> Vector{Any}

The genealogy path from `quantity` up to its root potential, as the list
of quantity *types* `[typeof(quantity), parent, …, FreeEnergy]`.  A root
(or non-genealogy) quantity returns the singleton `[typeof(quantity)]`.

```julia
differentiation_chain(SusceptibilityZZ())
# [SusceptibilityZZ, MagnetizationZ, FreeEnergy]  — i.e. χ ⟵ M ⟵ F
```

The chain terminates because the genealogy is a finite tree rooted at
the free energy; a cycle in the declared edges would loop forever and is
guarded against with an explicit depth cap.
"""
function differentiation_chain(q::AbstractQuantity)
    chain = Any[typeof(q)]
    cur = typeof(q)
    for _ in 1:32   # depth cap: the real tree is ≤ 3 deep; guards a mis-declared cycle
        e = derivative_edge(cur)
        e === nothing && return chain
        push!(chain, e.parent)
        cur = e.parent
    end
    return error("differentiation_chain: genealogy exceeded depth cap — cyclic edge?")
end
export differentiation_chain

"""
    potential_root(quantity) -> Type

The root potential of `quantity`'s genealogy — the last entry of its
[`differentiation_chain`](@ref).  For every response function in the
thermodynamic tree this is [`FreeEnergy`](@ref); everything ultimately
derives from the free energy (hence, via `F = −β⁻¹ ln Z`, from the
partition function).
"""
potential_root(q::AbstractQuantity) = last(differentiation_chain(q))
export potential_root

"""
    derivative_order(quantity, field::AbstractField) -> Int

How many times `quantity` is differentiated with respect to `field` on
the way up to the root — the order of `quantity` as a `field`-derivative
of its root potential.

```julia
derivative_order(SusceptibilityZZ(), MagneticField())  # 2  (χ = ∂²F/∂h²)
derivative_order(MagnetizationZ(),   MagneticField())  # 1  (M = ∂F/∂h)
derivative_order(SpecificHeat(),     MagneticField())  # 0  (no field derivatives)
```
"""
function derivative_order(q::AbstractQuantity, field::AbstractField)
    ft = typeof(field)
    n = 0
    cur = typeof(q)
    for _ in 1:32
        e = derivative_edge(cur)
        e === nothing && return n
        e.field === ft && (n += 1)
        cur = e.parent
    end
    return error("derivative_order: genealogy exceeded depth cap — cyclic edge?")
end
export derivative_order
