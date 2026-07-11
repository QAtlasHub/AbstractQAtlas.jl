"""
    AbstractQAtlas

The model-independent layer of the QAtlas ecosystem, in the spirit of
`AbstractFFTs`: concrete atlases (QAtlas) *implement* this package,
never the reverse.

Two responsibilities:

1. **core/** — the abstract type vocabulary for physical quantities:
   `AbstractQAtlasModel`, `BoundaryCondition` (`Infinite`/`OBC`/`PBC`),
   `AbstractQuantity` and its hierarchy, the generic
   [`fetch`](@ref AbstractQAtlas.fetch) verb, and the
   [`Universality`](@ref) machinery — so atlases and third packages
   share dispatch types without depending on a full atlas.

2. **relations/** — generic, model-independent physics relations as
   first-class tested objects with a uniform three-verb interface
   ([`residual`](@ref) / [`check`](@ref) / [`solve`](@ref)):
   critical-exponent scaling laws, fluctuation–dissipation identities,
   Wick's theorem, standard topological invariants, and finite-size
   scaling forms.

**Values do not live here.**  Reference numbers (critical temperatures,
exact magnetizations, exponent tables, …) belong to the implementing
atlas; this package owns only what is true independently of any model.
"""
module AbstractQAtlas

using LinearAlgebra: det, eigen, Hermitian

# core — type vocabulary
include("core/types.jl")
include("core/quantities.jl")
include("core/universality.jl")
include("core/distributions.jl")

# structure — model-independent definitional correspondences between the
# core quantities: transition classification and the critical
# quantity⇄exponent map from which the scaling forms are derived.
include("structure/transitions.jl")
include("structure/criticality.jl")

# relations — model-independent identities and forms
include("relations/interface.jl")
include("relations/scaling.jl")
include("relations/thermodynamic.jl")
include("relations/fundamental.jl")
include("relations/statistics.jl")
include("relations/wick.jl")
include("relations/topology.jl")

end # module AbstractQAtlas
