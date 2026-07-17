# Design: the type-keyed relation interface

Status: **Phase 1 landed** — the interface + the undecorated Green's-function
domain (Dyson, SpectralFromGreens, the six Keldysh relations) are migrated and
green; the remaining domains are pending (tracking issue #77). The end-state code
sketches below (e.g. §2's "`quantity_links.jl` is deleted") describe the *target*,
not the post-Phase-1 tree — deletion is **incremental**, one domain at a time, so
`quantity_links.jl` still holds the not-yet-migrated relations' entries.
This document is the source of truth for the migration; it supersedes the
scattered notes on #77 and folds in the review findings below.

## 1. Problem — symbol-keyed variables are a *correctness* hazard, not just an ergonomic one

Every relation is declared with formula-letter **symbols**:

```julia
@relation :spectral Dyson(G, G0, Σ)             = ...   # G  = RetardedGreensFunction
@relation :spectral SpectralFromGreens(A, ImGR) = ...   # ImGR from the same G
@relation :transport KelvinRelation(Π, S, T)    = Π - T*S   # S = thermopower (Seebeck)
@relation :fundamental FreeEnergyLegendre(F, U, S, β) = ... # S = thermal entropy
```

The symbol is chosen **independently of the quantity type**, so it can *drift*
(`RetardedGreensFunction` is `G` in `Dyson` but `GR` in `SpectralFromGreens`/
`Keldysh`) and it can *collide* (`:S` is used for **four** distinct quantities:
thermal entropy, thermopower, dynamical structure factor, and von Neumann
entropy; `:β` is used for both the inverse temperature and the critical
exponent).

This is not cosmetic. The registry, the `derive`/`derivable` solver, and the
derivation graph all key on the bare `Symbol`, so **a symbol collision fuses
physically-distinct quantities into one node and lets the verify-engine chain
across them**. Demonstrated on the current `main`:

```julia
julia> using AbstractQAtlas
# Pass S meaning THERMAL ENTROPY (3.0) and a temperature (2.0).
# There is no thermoelectric measurement anywhere in this data.
julia> AbstractQAtlas.derivable(S = 3.0, T = 2.0)
Set([:S, :T, :dF_dT, :dlnσ_dε, :spectral_integral, :Π])
#                    ^^^^^^^^  ^^^^^^^^^^^^^^^^^  ^^  all WRONG: S reused as
#                    thermopower / structure-factor, not entropy

julia> AbstractQAtlas.derive(:Π; S = 3.0, T = 2.0, debug = true)
DerivationTrace(:Π = 6.0  [indirect])
  1. KelvinRelation: {S, T} → :Π          #  Π = T·S, entropy silently used as Seebeck
```

The engine reports the Peltier coefficient `Π = 6.0` derived from a thermal
entropy — a dimensionally-meaningless number presented as a computed physical
quantity. This is exactly the "looks-like-success, wrong mechanism" failure the
verify-engine exists to prevent, occurring *inside* the verify-engine.

`quantities(rel)` — the relation → quantity map that powers
`relations_constraining` and the unified graph — is maintained **by hand** in
`relations/quantity_links.jl` (68 entries). It is a second, parallel encoding of
information the typed declaration would already carry, and therefore a standing
drift source.

**"One canonical symbol per type" cannot fix this**, because a type legitimately
appears under many symbols — the symbol carries a *decoration* the bare type does
not:

| type | symbols in use | decoration |
|---|---|---|
| `VonNeumannEntropy` | `S`, `S_A`, `S_B`, `S_AB`, `S_BC`, `S_ABC` | **subsystem / region** |
| `DynamicalStructureFactor` | `S`, `S_plus`, `S_minus` | **±ω** |
| `Resistivity` | `ρxx`, `ρxy` | **tensor component** |
| `SpecificHeat` | `C`, `Cp`, `Cv` | constraint (p / v) |
| `RetardedGreensFunction` | `G`, `GR` | **none — this one is the pure bug** |

The identity we actually want to key on is `(type, decoration)`, and the type
layer already exists (`AbstractQuantity`, `AbstractField`).

## 2. Target model — the *type* is the identity; formula letters are private

```julia
abstract type RelationVariable end
#   AbstractQuantity    observables      (RetardedGreensFunction, VonNeumannEntropy, Susceptibility{I}) — EXISTS
#   AbstractField       control fields   (InverseTemperature, Temperature, MagneticField, ChemicalPotential) — EXISTS
#   AbstractCoordinate  eval coordinates (Frequency, Momentum) — NEW, small
#   AbstractExponent    critical exps    (α, β, γ, ν, η, δ, z) — NEW, small

@relation :spectral SpectralFromGreens(A::SpectralFunction, G::RetardedGreensFunction) =
    A + imag(G) / π
# `A`, `G` are PRIVATE local bindings for the math; the exposed key is the TYPE.

bag = (RetardedGreensFunction => gr, SpectralFunction => a, SelfEnergy => σ)
relation_report(bag)      # matched by type, never by a string
```

`quantities(rel)` is **auto-derived** from the typed declaration, so
`quantity_links.jl` is **deleted**. `G`-vs-`GR`, the `:S` collision, and the
exponent-`β`-vs-inverse-temperature-`β` collision all become *structurally
impossible* (distinct types ⇒ distinct keys). `domain=` survives only for
*scoping* a report, never for *disambiguation*.

## 3. Capability-preservation map (break the keying, keep the physics)

| capability | how it is preserved |
|---|---|
| `residual` / `check` / `solve`, exact-arithmetic (`Rational` in ⇒ `Rational` out) | **residual math kernels stay verbatim**; only the key→value expansion wrapper changes ⇒ risk is confined to the keying layer, correctness held by existing tests |
| verify-engine (bag → every applicable identity) | type-keyed bag, identical semantics |
| derivation graph / solver (`derive`, `derivable`) | nodes `Symbol → RelationVariable type`; this graph and the quantity graph then **share a node space and unify into one typed network**. Traversal preserved and *enriched* (must-keep). |
| AD `thermal_derivative` | already dispatches on quantity types → nearly unchanged |
| inequalities (slack / `≥ 0` check) | unchanged |
| β-or-T | an `InverseTemperature ↔ Temperature` typed conversion, replacing the `Symbol`-level special case in `_normalize_kwargs` |

## 4. Review refinements (pressure-testing the approved shape)

These four points were surfaced reviewing the plan before committing to the
multi-PR migration. They are folded into §6.

**R1 — key on `(type, support)` from day one, never on bare type.**
The entropy relations are intrinsically multi-instance-of-one-type:
`I(A:B) = S(A) + S(B) − S(A∪B)` mentions `VonNeumannEntropy` three times with
three different regions. Bare-type keys *cannot* express this; only
`(type, support)` (or, equivalently, a decorated parametric type
`VonNeumannEntropy{RegionAB}`) can. If the interface keys on bare type for the
easy domains and adds `support` later, the entanglement/correlation domains force
a **second re-key of the whole interface**. So the variable key is
`(type, support)` *from the start*, with `support = global` as the default that
makes the undecorated domains read exactly like bare-type keys. This costs
nothing for spectral/transport and avoids double-churn.

**R2 — sequence the two changes; do not let the hard nut block the win.**
Type-keying (kill collisions, auto-derive `quantities`, unify the graph) is a
*preservation refactor* — its correctness is proved by the existing physics
tests going green under new keys. The `support` / `Region` layer
(entanglement auto-discovery, SSD `β(x)`, area law, lattice seam) is a genuinely
*new capability* with open research design (quantified relations over regions).
Land **all** of type-keying first, on the undecorated domains, delivering the
correctness fix. Then build `support`/`Region` as its own epic. Concretely,
`quantity_links.jl` is deleted **incrementally**: an entry is removed when its
domain migrates, and the entanglement entries stay until the region layer lands —
never delete it wholesale before the entropy domain can be keyed.

**R3 — make the identity/supplied split explicit; not everything is a type.**
Two kinds of variable must be distinguished in the declaration:

- **identity-bearing** (`AbstractQuantity` / `AbstractField` / `AbstractExponent`)
  — keyed by type, participate in matching, the graph, and `derive`.
- **supplied slots** (evaluation coordinates `ω`, `q`; supplied integrals
  `pv_imag`, `∫A` — the #50 functional boundary) — passed positionally/by name,
  *not* promoted to graph nodes. Coordinates recur (`ω₁, ω₂` in nonlinear
  response) and would need their own decoration, but they are evaluation points,
  not subjects of the identity, so they stay lightweight supplied parameters.

**R4 — the unified graph's node type broadens; keep a quantity-only view.**
After unification the node type is the broad `RelationVariable` union, so
`quantity_path` / `related_quantities` now traverse through field and exponent
nodes too. That enrichment is wanted, but the existing quantity-only queries
must keep a *projected view* (filter nodes to `AbstractQuantity`) so their
semantics do not silently change.

Not fixed by type-keying, tracked separately (see §7): the derivation graph's
`graph_reachable` **over-approximates** computability because it projects each
`{inputs} → output` hyperedge to one edge per input (AND collapsed to OR). That
is orthogonal to the key type and remains guarded by the `!!! warning` on
`derivation_graph` and by keeping `derivable`/`derive` (which require *all*
inputs and call the real `solve`) as the honest path.

## 5. `support`, `Region`, and the lattice seam (the second epic)

A variable is `(type, support, extra)`. `support` = *where / on what it is
evaluated*:

| support | example |
|---|---|
| `global` (default) | a bulk observable |
| **Region** (a set of sites) | `VonNeumannEntropy(region)` |
| point `x` / profile `f(x)` | **SSD** `InverseTemperature(at = x)` / `β(x)` (modular Hamiltonian, sine-square deformation) |
| pair `(r, r′)` | correlation `⟨O(r)O(0)⟩` |

Entropy regions and SSD position-dependent fields collapse into one abstraction.
The **set layer** of `support` carries the entropy inequalities; the **geometric
layer** carries the area law.

`Region = Set{AbstractSite}` with `∪` / `∩` / disjointness, written over abstract
regions. Payoff — the verify-engine gets *smarter*: drop entropies of many
concrete regions into a bag and the engine **auto-discovers** every
subadditivity / SSA / Araki–Lieb instance whose region unions are present, with
no `A` / `B` / `AB` hand-labeling.

**Dimension-agnostic + optional lattice.** `AbstractSite` is abstract
(1D `Int` / 2D `(Int,Int)` / ND `NTuple{D,Int}` / a label), so the set algebra
and every entropy inequality are **ND from day one, stdlib-only**. Dimension
enters *only* the geometric layer (contiguity, boundary `∂A`, distance, SSD
profile `sin²(πx/L)`), which is an **optional lattice extension** (inject a
lattice; `LatticeCore` is one provider, a 1D chain or a future ND lattice are
others). The core stays lattice-free and ND-ready regardless of `LatticeCore`'s
2D focus.

## 6. Migration plan (functionality-preserving, multi-PR)

1. **Do not touch the residual math kernels.** Physics correctness is held by
   the existing tests; risk is confined to keying.
2. Define the `RelationVariable` layer: reuse `AbstractQuantity` /
   `AbstractField`; add `AbstractCoordinate`, `AbstractExponent`. Design the
   variable key as `(type, support)` with `support = global` default (R1).
3. Rebuild `@relation` / `@inequality` to take **typed** variables, keeping the
   residual kernels verbatim and auto-generating `quantities(rel)`. Keep the
   identity/supplied split (R3).
4. Rebuild `residual` / `check` / `solve` over type keys (preserve exact
   arithmetic + the generic affine `solve`). β-or-T becomes the
   `InverseTemperature ↔ Temperature` typed conversion.
5. Rebuild `relation_report` / `applicable_relations` / `check_all` as a
   **type-keyed bag**.
6. **Prototype end-to-end on the spectral / Green's domain** (undecorated
   singletons, no regions): validate keying, the unified network traversal, and
   the verify-engine before rolling out.
7. Roll out **per domain, big-bang per domain**, on the undecorated domains
   first (spectral, transport, thermodynamic, scaling, quantum, wick, ensemble).
   Delete each domain's `quantity_links.jl` entries as it migrates (R2).
8. Unify the derivation graph (now typed nodes) with the quantity graph into one
   `KnowledgeGraph`; preserve `derive` / `derivable` / path-finding and add the
   quantity-only projected view (R4). Network traversal is a must-keep.
9. Migrate tests: **physics assertions unchanged, only the keying syntax
   changes**. Full green is the preservation proof.

Then, as a separate epic (§5):

10. Build `support` / `Region` (dimension-agnostic, stdlib-only) and migrate the
    entanglement + correlation-pair domains onto `(type, Region)` /
    `(type, pair)`; delete the remaining `quantity_links.jl` entries.
11. Add the geometry seam (boundary / area-law / SSD profile) via an optional
    lattice extension; keep the core lattice-free.

## 7. Open questions to settle during design

- **Role/decoration encoding**: parametric type `VonNeumannEntropy{Region}`
  (dispatchable, consistent with `Susceptibility{I}`) vs a `(type, tag)` pair.
  Leaning parametric — it makes distinct regions distinct *types*, so the
  type-keying machinery covers the multi-instance case directly.
- **How far to type non-observables**: exponents & fields → typed (kills the `β`
  collision); coordinates (`ω`, `q`) & supplied integrals (`pv_imag`, `∫A`) →
  lightweight explicit "supplied" slots, not graph nodes (R3).
- **`Region`**: atom granularity (sites vs named blocks), disjointness detection,
  normalization for matching in auto-discovery.
- **Derivation-graph over-approximation** (§4, orthogonal): whether to make
  reachability hyperedge-aware (require all inputs) or keep the current
  structural-graph + honest-`derivable` split. Not a blocker.

## 8. Parametric tensor quantities — the next-wave resolution

Transport, and the tensor response quantities generally, are index-parametric —
`Conductivity{I}`, `Susceptibility{I}`, `Thermopower{I}`, `Resistivity{I}`,
`ThermalConductivity{I}`, `DrudeWeight{I}`, `Magnetization{A}`,
`DynamicalSusceptibility{I}`, … The Phase-1 review flagged (finding C4) that the
bag matches by `===` while the quantity graph erases to `_family`, so a slot keyed
on the bare family `Conductivity` can never `===`-match a concrete
`Conductivity{(:x,:x)}` bag entry — and the load-time guard now rejects such a
bare-family declaration loudly. This section resolves how the parametric domains
migrate.

**Three index structures appear across the relations:**

| structure | example | slots |
|---|---|---|
| same index on every slot | `WiedemannFranz` κ_ii = L₀ σ_ii T | `ThermalConductivity{ii}`, `Conductivity{ii}` |
| fixed, *different* indices | `HallAngle` tanθ = σ_xy / σ_xx | `Conductivity{(:x,:y)}`, `Conductivity{(:x,:x)}` |
| index *unification* across slots | `SusceptibilityFDT` χ_AB = β Cov(M_A, M_B) | `Susceptibility{(A,B)}`, `Magnetization{A}`, `Magnetization{B}` |

**Key insight — concrete component keys work with the SHIPPED machinery, today.**
A concrete instantiation like `Conductivity{(:x,:y)}` is `isconcretetype`, is
`<: RelationVariable`, and `Conductivity{(:x,:y)} !== Conductivity{(:x,:x)}` — so
distinct components are distinct `VariableKey`s that pass both load guards
(concrete ✓, distinct ✓) and match by `===` exactly. So the **fixed-index** and
**same-index-when-written-concretely** cases (HallAngle, the Hall/resistivity
family, and any per-component check) migrate **now, with no new machinery** — each
relation is simply component-specific, which is physically what HallAngle *is*.

**Phased plan:**

- **Phase 1 (no new machinery):** migrate the parametric domains declaring
  **concrete component types** where the component is fixed. Component-agnostic
  scalar laws (WF, Mott, …) are declared on their canonical (longitudinal)
  component. `_auto_quantities` must `unique`-dedup (a relation with two
  `Conductivity{…}` slots erases to `(Conductivity, Conductivity)`; the old
  hand-link was `(Conductivity,)`).
- **Phase 2 (the real §7 — deferred):** **family-generic slots with index
  unification**, so the verify-engine *auto-discovers* every component instance
  from a bag (drop `χ_xy, M_x, M_y` → find `χ_xy = β Cov(M_x, M_y)` without a
  hand-written per-component relation). This needs (a) macro support for a type
  variable — `SusceptibilityFDT(χ::Susceptibility{(A,B)}, M_A::Magnetization{A},
  M_B::Magnetization{B}) where {A,B}` — and (b) a bag matcher that unifies the
  variable across slots against the concrete keys present. This is the tensor
  analogue of the entanglement `Region` auto-discovery (§5–§6) and shares its
  "quantified relations" machinery; build them together.

The load guard is the safety net between the phases: until Phase 2 lands, a
bare-family declaration fails at load rather than silently matching nothing.

## 9. Status of the in-flight symbol-based work

PR #78 (the symbol-based D1: complex `SpectralFromGreens`, `Dyson` `:G → :GR`
alignment, bag-adoption test) and the `feat/spectral-bag-turnkey` branch are
**superseded / paused** by this redesign — the `G` / `GR` fix is subsumed (both
become `RetardedGreensFunction`). Do not invest further in the symbol-keyed
approach.
