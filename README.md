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

| lives here | lives in the implementing atlas |
|---|---|
| type vocabulary — `AbstractQAtlasModel`, `AbstractQuantity` (+ hierarchy), `BoundaryCondition` (`Infinite`/`OBC`/`PBC`), the generic `fetch` verb, `Universality{C}` | concrete models and registered `fetch` methods |
| generic **relations** — scaling laws, fluctuation–dissipation identities, Wick's theorem, topological invariants, FSS forms | reference **values** — critical temperatures, exact magnetizations, exponent tables |

A relation is an *identity among observables or exponents* — true
independently of any model.  Expressing each one once, as a tested
object, means downstream packages stop re-deriving them ad hoc in
comments and per-model tests.

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
- **FSS forms**: `collapse_coordinates`, `fss_peak_scaling`,
  `order_parameter_form`, `correlation_length_form`,
  `susceptibility_form`

Every relation is tested against an *independent* expectation: exact
rational exponent sets, derivative-vs-fluctuation cross-checks,
Fock-space ED vs the Wick determinant, and known topological phase
structures (SSH, QWZ).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sotashimozono/AbstractQAtlas.jl")
```

## Consumers

- [QAtlas.jl](https://github.com/sotashimozono/QAtlas.jl) — declares the
  reference values; adopts this package as its type + relations base.
- [ClassicalMonteCarlo.jl](https://github.com/sotashimozono/ClassicalMonteCarlo.jl)
  — finite-size-scaling validation consumes the FSS forms and gates
  extracted exponents with the scaling relations.

## Roadmap

- Full migration of the remaining QAtlas quantity tags (NMR, Loschmidt,
  entanglement family, structure factors, …).
- Anomalous / BdG Wick contraction (Pfaffian).
- More identities (Maxwell relations, Kramers–Kronig for response
  functions).
