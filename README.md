# AbstractQAtlas.jl

[![docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://qatlashub.github.io/AbstractQAtlas.jl/stable/)
[![docs: dev](https://img.shields.io/badge/docs-dev-purple.svg)](https://qatlashub.github.io/AbstractQAtlas.jl/dev/)
[![Build Status](https://github.com/QAtlasHub/AbstractQAtlas.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/QAtlasHub/AbstractQAtlas.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/main/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

*The model-independent core of the QAtlas ecosystem.*

**QAtlas** is a knowledge base of quantum many-body physics: for a given model —
the transverse-field Ising chain, the Hubbard model, and so on — it records
reference values of physical quantities. **AbstractQAtlas** is the layer beneath
it, holding the part of that knowledge that does *not* depend on the model. The
numbers live in QAtlas; the vocabulary of quantities and the universal laws
relating them live here.

Three ideas run through the package:

- **Quantity** — a physical observable, written as a type: `Energy`,
  `Susceptibility`, `SpectralFunction`, …
- **Domain** — the area of physics a quantity belongs to: thermodynamics,
  criticality, correlations, transport, quantum information, topology. Every
  quantity is filed under one.
- **Relation** — a universal identity connecting quantities, true for *any*
  model: `χ = ∂M/∂h`, the fluctuation–dissipation theorem, the critical scaling
  laws. Each is a first-class, tested object.

In the spirit of `AbstractFFTs`, the concrete atlas
([QAtlas.jl](https://github.com/QAtlasHub/QAtlas.jl)) *implements* this package;
the dependency never points the other way.

## Quantities, by domain

| domain | quantities |
|---|---|
| **thermodynamics** | `Energy`, `FreeEnergy`, `ThermalEntropy`, `SpecificHeat`, `PartitionFunction`, `Pressure` |
| **criticality** | critical exponents `(α, β, γ, δ, ν, η)`, `CorrelationLength`, `MassGap` |
| **correlations** | `RetardedGreensFunction`, `KeldyshGreensFunction`, `SpectralFunction`, `DynamicalStructureFactor`, `SelfEnergy` |
| **transport** | `Susceptibility` (linear *and* nonlinear), `Conductivity`, `ThermalConductivity`, `Thermopower`, `HallCoefficient` |
| **quantum information** | `VonNeumannEntropy`, `RenyiEntropy`, `MutualInformation`, `Concurrence`, `TopologicalEntanglementEntropy` |
| **quantum foundations** | `KineticEnergy`, `PotentialEnergy`, `EnergyVariance` |
| **topology** | `ChernNumber`, `BerryCurvature`, `BoundaryModeCount` |

Tensor quantities keep their indices: `Susceptibility(:x, :y)` is the
off-diagonal `χ_xy`, and `Susceptibility(:x, :y, :z)` the second-order
`χ⁽²⁾_{x;yz}` — not blurred into a scalar.

## Relations

Every relation carries three verbs:

- `residual` — how far the values are from satisfying it,
- `check` — whether they satisfy it,
- `solve` — rearrange it for one unknown.

```julia
using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

solve(Widom(), Val(:δ); β = 1//8, γ = 7//4)               # 15//1  — a scaling law
check(SusceptibilityFDT(); χ = 0.5, var_M = 2.0, T = 4.0)  # χ = β·Var(M)  ?
residual(Dyson(); G = g, G0 = g0, Σ = σ)                  # G⁻¹ = G₀⁻¹ − Σ
```

The catalogue spans every domain above — the thermodynamic web (`F = U − TS`,
the Maxwell relations), the scaling laws, the fluctuation–dissipation theorem
(thermodynamic, dynamical, and Keldysh `G^K = h(ω)(G^R − G^A)`), Kramers–Kronig,
Wick's theorem, the Kubo formula, and the topological invariants. Rational
inputs give rational residuals, so an exactly known value satisfies its relation
*exactly*, not just to floating-point tolerance.

## Following relations from one quantity to another

The relations connect quantities, so together they form a graph — and the graph
is queryable. Ask how two quantities relate, or hand over what you measured and
let the package derive another by finding a route:

```julia
using AbstractQAtlas: related_quantities, quantity_path, derive

related_quantities(Susceptibility(:z, :z))
#   Magnetization          — as ∂/∂h, and again through fluctuation–dissipation
#   StaticStructureFactor  — through the structure-factor sum rule

quantity_path(SpecificHeat(), Magnetization(:z))
#   SpecificHeat — Energy — FreeEnergy — Magnetization      (their common root)

derive(:δ; α = 0//1, β = 1//8, debug = true)
#   DerivationTrace(:δ = 15//1  [indirect])
#     1. Rushbrooke: {α, β} → :γ
#     2. Widom:      {β, γ} → :δ
```

`derive` chains equalities only, calls the real `solve` at each step, and with
`debug = true` reports the exact route — so an indirectly obtained value is
auditable, not taken on trust.

## Citations, managed

Every relation drawn from the literature keeps its source. The references sit in
`docs/references.bib`, and CI checks each DOI against Crossref: a citation is
either real and resolvable or absent — never invented. Kubo (1957) for the
response functions, Callen–Welton (1951) for the FDT, Sugiura–Shimizu for the
thermal-pure-quantum construction, and so on.

## Installation

```julia
using Pkg
Pkg.add("AbstractQAtlas")
```

Depends only on Julia standard libraries. `ForwardDiff` is a weak dependency;
loading it turns on automatic differentiation of the response genealogy.

## Ecosystem

- [QAtlas.jl](https://github.com/QAtlasHub/QAtlas.jl) — the concrete atlas that
  supplies the reference values and implements this package's types and
  relations.

Beyond backing an atlas, AbstractQAtlas is meant as a verification engine: a
calculation (ED, MPS, TPQ, …) can check its measured quantities straight against
the relation web — fluctuation–dissipation, ensemble equivalence, scaling, Kubo
— instead of re-deriving each identity by hand.
