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

**Values do not live here, and neither do model-specific laws.**  This
package owns only what holds **universally within a domain** — independent
of the system's symmetry, Hamiltonian, or any individual detail (Wick's
theorem, the fluctuation–dissipation and Maxwell relations, the entropy
inequalities, Kramers–Kronig, the scaling laws, …).  Reference *numbers*
(critical temperatures, exact magnetizations, exponent tables) AND
*model-specific* relations (the Drude mobility `μ=eτ/m`, the ±J Nishimori-
line energy, the SK de Almeida–Thouless line, single-band `R_H=1/ne`, …)
belong to the implementing atlas (QAtlas), not here.  The library is a
**universal yardstick**: apply its relations to measured quantities to
check whether a system obeys the laws that must hold regardless of its
details.
"""
module AbstractQAtlas

using LinearAlgebra: det, eigen, Hermitian
using InteractiveUtils: subtypes    # quantity-graph node discovery (structure/graph.jl)

# core — type vocabulary
include("core/types.jl")
include("core/indices.jl")
include("core/representations.jl")
include("core/quantities.jl")
include("core/universality.jl")
include("core/distributions.jl")
include("core/fields.jl")
include("core/relation_variables.jl")   # the RelationVariable layer (type-keyed variables)
include("core/region.jl")                # the Region set layer (entanglement support, §5)

# structure — model-independent definitional correspondences between the
# core quantities: transition classification, the critical
# quantity⇄exponent map from which the scaling forms are derived, and the
# response-function derivative genealogy rooted at the free energy.
include("structure/transitions.jl")
include("structure/criticality.jl")
include("structure/potentials.jl")        # thermodynamic potentials → the Maxwell relations
include("structure/scaling_dimensions.jl")
include("structure/response.jl")
include("structure/tensor_symmetry.jl")
include("structure/spectral.jl")
include("structure/keldysh.jl")
include("structure/fourier.jl")

# relations — the CORE relation machinery (the AbstractFFTs-like interface:
# the @relation / @inequality macros, residual/check/solve/slack, the
# registry, relation_report/applicable_relations) lives in the parent…
include("relations/interface.jl")

# …and the physics is organized into submodules by PHYSICAL DOMAIN (not by
# file): each `using ..AbstractQAtlas` for the macro + shared vocabulary,
# includes its (bare) relation files, and is flat-re-exported below — so
# `using AbstractQAtlas` and `AbstractQAtlas.Name` keep working while the
# source reflects the physics.  The fine-grained per-relation `domain` tag
# is preserved (for `relation_report` filtering) INSIDE each module.

"Equilibrium statistical mechanics: occupation statistics → ensembles → thermodynamic potentials → response, FDT & stability."
module StatisticalMechanics
    using ..AbstractQAtlas
    using ..AbstractQAtlas: _beta          # β-or-T normalization (occupation functions)
    import ..AbstractQAtlas: _solve        # extended for a non-affine variable (FreeEnergyFromZ:Z)
    include("relations/thermodynamic.jl")
    include("relations/ensembles.jl")
    include("relations/fluctuation.jl")
    include("relations/statistics.jl")
    include("relations/fundamental.jl")
end

"Critical phenomena and conformal field theory: scaling laws, finite-size scaling, Cardy, the c-theorem."
module Criticality
    using ..AbstractQAtlas
    include("relations/scaling.jl")
    include("relations/cft.jl")
end

"Correlations, Green's functions and response: the spectral graph (Dyson, A=−ImG/π), the Keldysh RAK structure + fluctuation–dissipation, Wick / Bloch–De Dominicis (Gaussian factorization), Kramers–Kronig, detailed balance."
module Correlations
    using ..AbstractQAtlas
    using LinearAlgebra: inv, det
    include("relations/spectral.jl")
    include("relations/keldysh.jl")
    include("relations/wick.jl")
end

"Transport: DC/AC conductivity, thermal & thermoelectric coefficients, the Hall family, Onsager, Wiedemann–Franz, optical sum rule, Johnson–Nyquist."
module Transport
    using ..AbstractQAtlas
    include("relations/transport.jl")
end

"Quantum information & entanglement: the entropy zoo, its inequalities, multipartite entanglement, measurement and topological entanglement entropy."
module QuantumInformation
    using ..AbstractQAtlas
    include("relations/entanglement.jl")
end

"Quantum-mechanical foundations & bounds: virial, Hellmann–Feynman, Ehrenfest, zero-variance eigenstate, the uncertainty relation and the Lieb–Robinson bound."
module QuantumFoundations
    using ..AbstractQAtlas
    include("relations/quantum.jl")
end

"Topological invariants: Chern number, TKNN, winding, bulk–boundary correspondence."
module Topology
    using ..AbstractQAtlas
    using LinearAlgebra: det, eigen, Hermitian
    include("relations/topology.jl")
end

# flat re-export: lift every physics submodule's public names to the top
# level (stdlib-only — no Reexport.jl), preserving the flat API.
const _PHYSICS_MODULES = (
    StatisticalMechanics,
    Criticality,
    Correlations,
    Transport,
    QuantumInformation,
    QuantumFoundations,
    Topology,
)
for _M in _PHYSICS_MODULES
    for _n in names(_M)
        _n === nameof(_M) && continue
        @eval using .$(nameof(_M)): $_n
        @eval export $_n
    end
end

# the relation → quantity map (needs every relation + quantity name at the
# top level, so it goes AFTER the re-export) — makes the registry queryable
# via `quantities` / `relations_constraining`.
include("relations/quantity_links.jl")

# region-keyed entanglement-entropy auto-discovery (needs the re-exported
# Subadditivity / ArakiLieb + the Region set layer): region_report / region_check_all.
include("relations/region_entropy.jl")

# the abstract typed-graph PARENT: one KnowledgeGraph{N} kernel (nodes, typed
# edges, traversal, reachability, shortest-path, JSONL export) that every
# concrete graph below — quantity, derivation, and QAtlas's model graph — is an
# instance of.  Generic, stdlib-only, no rendering.
include("graph/kernel.jl")

# the quantity-relationship graph as a KnowledgeGraph{Type}: folds the separate
# typed edges (derivative_edge / spectral_origin / fourier_conjugate + the
# relation↔quantity links) into one graph — related_quantities / quantity_path /
# quantity_graph_jsonl.  Needs every edge accessor AND the relation↔quantity map.
include("structure/graph.jl")

# the DIRECTED derivation graph + lazy solver: reads the registry as a graph
# whose edges are `solve` directions, so a target quantity can be derived from
# known ones by finding and running one route (`derive` / `derivable`), with a
# debug trace naming the indirect path for safety.  Needs the full registry.
include("relations/derivation.jl")

# the ecosystem report/card contract: `report(model, quantity, bc; value, …)` packages
# an oracle's computed value into a schema-v2 `Card` (the reporter-facing sibling of
# `fetch`), streamable via `card_jsonl`.  Reporters depend on AbstractQAtlas only and
# push DATA — cycle-free.  Needs the core type vocabulary + graph/kernel's `_json_str`.
include("report/card.jl")

# the functional-evaluation interface (the scope-line #14 seam): generic verbs
# `principal_value_hilbert` / `spectral_moment` over an `AbstractResponse`, owned here
# (like `fetch`), evaluated by the functional sibling — so the supplied-integral relations
# become turnkey from a measured response once that package is present.
include("evaluation.jl")

# automatic-differentiation entry point (methods live in ext/, ForwardDiff)
include("autodiff.jl")

end # module AbstractQAtlas
