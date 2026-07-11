# core/types.jl — the root type vocabulary of the QAtlas ecosystem.
#
# Ported from QAtlas `src/core/type.jl` MINUS the deprecation surface
# (`Model{M}` / `Quantity{Q}` phantom wrappers, `AbstractModel` alias,
# `canonicalize_*` shims) — legacy compatibility stays in QAtlas'
# `src/deprecate/`.  This file owns only the clean, forward-looking
# contract that concrete atlases implement.

"""
    AbstractQAtlasModel

Abstract parent type for every atlas model.  Concrete subtypes carry
their physics parameters as typed fields, e.g.

```julia
struct TFIM <: AbstractQAtlasModel
    J::Float64
    h::Float64
end
```

Implementing packages (QAtlas and friends) subtype this and register
[`fetch`](@ref) methods per `(model, quantity, bc)` triple.
"""
abstract type AbstractQAtlasModel end
export AbstractQAtlasModel

"""
    BoundaryCondition

Abstract parent type.  The three concrete subtypes carry system-size
information where applicable, so [`fetch`](@ref) can read it from the
BC instead of `kwargs`:

- [`Infinite`](@ref) — thermodynamic limit; no size.
- [`PBC`](@ref)`(N::Int)` — periodic boundary conditions at finite `N`.
- [`OBC`](@ref)`(N::Int)` — open boundary conditions at finite `N`.

For backward compatibility, the zero-argument constructors `PBC()` and
`OBC()` exist and set `N = 0`, which signals "caller will pass `N` via
kwargs" — legacy fetch methods still look at `kwargs[:N]`.  New fetch
methods read `bc.N` directly.
"""
abstract type BoundaryCondition end
export BoundaryCondition

"""
    Infinite()

Thermodynamic-limit boundary condition — no finite size.
"""
struct Infinite <: BoundaryCondition end
export Infinite

"""
    OBC(N::Int)
    OBC(; N::Int = 0)

Open boundary condition.  `N` is the chain length.  `N = 0` is a legacy
sentinel meaning "size unspecified — caller passes it via kwargs";
`fetch` methods that accept `OBC(0)` must look up `kwargs[:N]`.
"""
struct OBC <: BoundaryCondition
    N::Int
end
OBC(; N::Int=0) = OBC(N)
export OBC

"""
    PBC(N::Int)
    PBC(; N::Int = 0)

Periodic boundary condition.  See [`OBC`](@ref) for the `N = 0`
sentinel.
"""
struct PBC <: BoundaryCondition
    N::Int
end
PBC(; N::Int=0) = PBC(N)
export PBC

"""
    _bc_size(bc::BoundaryCondition, kwargs) -> Int

Return the effective system size for `bc`.  Prefers `bc.N` when it is
positive; otherwise looks up `kwargs[:N]`; otherwise throws.  Legacy
fetch methods can use this helper to accept both `OBC(N=24)` and
`OBC(); N=24` call forms.
"""
function _bc_size end

function _bc_size(::Infinite, kwargs)
    return error("_bc_size called on Infinite; call a size-free fetch method instead")
end
function _bc_size(bc::Union{OBC,PBC}, kwargs)
    bc.N > 0 && return bc.N
    haskey(kwargs, :N) && return Int(kwargs[:N])
    return error("$(typeof(bc)): N unspecified — pass via OBC(N=...) or kwargs N=...")
end

"""
    AbstractQuantity

Abstract parent type for quantities.  Concrete quantity structs
(e.g. `struct MagnetizationX <: AbstractQuantity end`) make dispatch
static and naming explicit (axis, entropy variant, …).
"""
abstract type AbstractQuantity end
export AbstractQuantity

"""
    fetch(model, quantity, bc; kwargs...)

Return the stored / computed value of `quantity` for `model` under
boundary condition `bc`.  The canonical signature takes a concrete
model struct + concrete quantity struct + BC.

`AbstractQAtlas` owns the *generic function* only; each implementing
package registers one method per supported `(model, quantity, bc)`
triple.  This top-level fallback throws an informative error for
un-implemented triples.
"""
function fetch(
    model::AbstractQAtlasModel, quantity::AbstractQuantity, bc::BoundaryCondition; kwargs...
)
    return error(
        "no fetch method for model=$(typeof(model)), " *
        "quantity=$(typeof(quantity)), bc=$(typeof(bc)). " *
        "The implementing package must define " *
        "`fetch(::$(typeof(model)), ::$(typeof(quantity)), ::$(typeof(bc)); ...)`.",
    )
end
# NOTE: `fetch` is deliberately NOT exported — it would clash with
# `Base.fetch`.  Call it qualified (`AbstractQAtlas.fetch`) or import it
# explicitly (`using AbstractQAtlas: fetch`), exactly as QAtlas does.
