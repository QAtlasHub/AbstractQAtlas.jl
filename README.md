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

## The three-verb interface

```julia
using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

# 2D Ising exponents are exact rationals — relations hold EXACTLY
# (Rational in ⇒ Rational out; no silent float promotion):
residual(Rushbrooke(); α=0//1, β=1//8, γ=7//4)   # 0//1
solve(Widom(), Val(:γ); β=1//8, δ=15//1)          # 7//4
check(Fisher(); γ=7//4, ν=1//1, η=1//4)           # true

# gate a whole exponent table at once:
exponents_consistent((α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4); d=2)

# fluctuation–dissipation: solve() doubles as the estimator downstream
# codes should use, so the formula lives in exactly one place:
#   c_v = β² Var(E) / N
solve(SpecificHeatFDT(), Val(:C); var_E=var_E, β=β, N=N)
```

Covered in v0.1:

- **Scaling laws**: `Rushbrooke`, `Widom`, `Fisher`, `Josephson`,
  `exponents_consistent`
- **Fundamental equations**: `FreeEnergyFromZ`, `FreeEnergyLegendre` (F = U − TS), `EntropyResponse`, `GibbsHelmholtz`
- **Thermodynamic identities**: `SpecificHeatFDT`, `SusceptibilityFDT`,
  `LinearResponseFDT`
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
