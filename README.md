# AbstractQAtlas.jl

[![docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://codes.sota-shimozono.com/AbstractQAtlas.jl/stable/)
[![docs: dev](https://img.shields.io/badge/docs-dev-purple.svg)](https://codes.sota-shimozono.com/AbstractQAtlas.jl/dev/)
[![Julia](https://img.shields.io/badge/julia-v1.11+-9558b2.svg)](https://julialang.org)
[![Code Style: Blue](https://img.shields.io/badge/Code%20Style-Blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

<a id="badge-top"></a>
[![codecov](https://codecov.io/gh/sotashimozono/AbstractQAtlas.jl/graph/badge.svg?token=Q3oEEiz9A2)](https://codecov.io/gh/sotashimozono/AbstractQAtlas.jl)
[![Build Status](https://github.com/sotashimozono/AbstractQAtlas.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sotashimozono/AbstractQAtlas.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/main/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**The model-independent layer of the QAtlas ecosystem** — abstract
quantity vocabulary + generic physics relations as first-class, tested
objects.  Zero non-stdlib dependencies.

In the spirit of `AbstractFFTs`: concrete atlases
([QAtlas.jl](https://github.com/sotashimozono/QAtlas.jl)) *implement*
this package, never the reverse.

Three layers, each model-independent:

- **`core/` — the nouns.** `AbstractQAtlasModel`, `AbstractQuantity` (+ hierarchy), `BoundaryCondition` (`Infinite`/`OBC`/`PBC`), the generic `fetch` verb, `Universality{C}`, the distribution vocabulary.
- **`structure/` — the definitions.** The generic *correspondences* between quantities that hold by definition: which critical exponent governs which observable (`Susceptibility ↦ γ`, `SpontaneousMagnetization ↦ β`), the phase-transition classification (first-order / continuous / BKT). From these, the singular and finite-size-scaling forms are *derived* — `χ_max ∼ L^{γ/ν}` falls out of `critical_scaling(Susceptibility)`, it is not a number you pass by hand.
- **`relations/` — the checks.** Identities among observables/exponents as first-class, tested objects (scaling laws, fluctuation–dissipation, Wick, topological invariants), verified against numbers via `residual`/`check`/`solve`.

The implementing atlas holds the **values** — critical temperatures, exact magnetizations, exponent tables. This package holds only what is true *independently of any model*: the vocabulary, the definitional correspondences, and the identities they must satisfy.

## Declare once, derive everything

A relation is written **exactly once** — one line, one mathematical
statement:

```julia
@relation :scaling Rushbrooke(α, β, γ) = α + 2β + γ - 2
```

That single declaration yields the struct, the residual kernel, the
`variables`/`domain` introspection traits, registry membership, and —
with **no hand-written rearrangements** — `solve` for every variable the
expression is affine in (three exact kernel probes; a non-affine
variable is *refused*, never silently mis-solved).  The β-or-T keyword
convention is normalized once, in the verb layer.

```julia
using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

# Rational in ⇒ Rational out; exactly-known values satisfy EXACTLY:
residual(Rushbrooke(); α=0//1, β=1//8, γ=7//4)   # 0//1
solve(Widom(), Val(:δ); β=1//8, γ=7//4)           # 15//1 — derived, not hand-coded
check(Fisher(); γ=7//4, ν=1//1, η=1//4)           # true

# fluctuation–dissipation: solve() doubles as the downstream estimator
solve(SpecificHeatFDT(), Val(:C); var_E=var_E, T=T, N=N)   # c_v = β²Var(E)/N
```

## Adopting from another package: one call

A consumer never hand-lists relations.  Hand over a `NamedTuple` of
whatever it measured or declares; every applicable relation is selected
by its variables and checked:

```julia
# gate an exponent table (atlas registry, MC-extracted exponents, …):
check_all((α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4, d=2))   # true, exactly

# cross-check measured thermodynamics, with per-relation diagnostics:
relation_report((C=c, var_E=v, T=T, N=N); atol=tol)

# discover what would be checked:
applicable_relations((F=f, U=u, S=s, β=β))    # → [FreeEnergyLegendre()]
```

An empty match is `false`, never a silent green.  Pass `domain=` when a
data set mixes families (physics overloads names: the exponent `β` vs
the inverse temperature `β`).  Downstream packages declare their own
relations with the same `@relation` macro — methods land on these
generics, and the verbs and report machinery apply unchanged.

Covered in v0.1:

- **Scaling laws**: `Rushbrooke`, `Widom`, `Fisher`, `Josephson`,
  `exponents_consistent`
- **Fundamental equations**: `FreeEnergyFromZ`, `FreeEnergyLegendre` (F = U − TS), `EntropyResponse`, `GibbsHelmholtz`
- **Thermodynamic identities**: `SpecificHeatFDT`, `SusceptibilityFDT`,
  `LinearResponseFDT`
- **Distributions & statistics**: `MicroCanonical` / `Canonical` / `GrandCanonical` / `Squeezed` + `ensemble_weight`; `Fermionic`/`Bosonic` occupation functions; `ThermalAverage` marker (⟨Q⟩_D at the type level)
- **Wick's theorem**: `wick_contraction`, `wick_density_correlation`
  (number-conserving Gaussian fermions)
- **Topological invariants**: `winding_number` (1D two-band),
  `chern_number` (Fukui–Hatsugai–Suzuki lattice method), `TKNN`

**Structure — definitional correspondences** (the `structure/` layer):

- **Critical correspondence**: `critical_scaling(quantity)` maps each
  observable to the exponent governing its singularity
  (`Susceptibility ↦ γ`, `SpontaneousMagnetization ↦ β`, `SpecificHeat
  ↦ α`, `CorrelationLength ↦ ν`; plus `critical_isotherm ↦ δ`,
  `correlation_decay ↦ η`). From it the forms are *derived*:
  `singular_form`, `fss_size_exponent` (`χ_max ∼ L^{γ/ν}`),
  `fss_peak`, `collapse_coordinates` — the exponent combination is
  looked up, never hand-passed.
- **Transition classification**: `FirstOrder` / `ContinuousTransition`
  / `KosterlitzThouless` with `ehrenfest_order`, `has_order_parameter`,
  `has_latent_heat`, `has_critical_exponents`.

Every relation is tested against an *independent* expectation: exact
rational exponent sets, derivative-vs-fluctuation cross-checks,
Fock-space ED vs the Wick determinant, known topological phase
structures (SSH, QWZ), and the correspondence reproducing the exact
textbook power laws.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sotashimozono/AbstractQAtlas.jl")
```

## Consumers

- [QAtlas.jl](https://github.com/sotashimozono/QAtlas.jl) — declares the
  reference values; adopts this package as its type + relations base.
- [ClassicalMonteCarlo.jl](https://github.com/sotashimozono/ClassicalMonteCarlo.jl)
  — finite-size-scaling validation consumes the critical correspondence
  and gates
  extracted exponents with the scaling relations.

## Roadmap

- Full migration of the remaining QAtlas quantity tags (NMR, Loschmidt,
  entanglement family, structure factors, …).
- Anomalous / BdG Wick contraction (Pfaffian).
- More identities (Maxwell relations, Kramers–Kronig for response
  functions).
