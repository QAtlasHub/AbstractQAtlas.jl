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

**Use as a verification engine.** Beyond serving an atlas, this is meant to be the relation web a quantum-many-body calculation (ED, MPS, TPQ, …) checks its *measured* quantities against — you compute an observable, then use the library's relations directly to test consistency (fluctuation–dissipation, ensemble equivalence, scaling laws, Kubo, …) rather than re-deriving each identity by hand. (Bridging these continuum relations to discrete grids / finite size is tracked in [#19](https://github.com/sotashimozono/AbstractQAtlas.jl/issues/19).)

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
- **Topological invariants & characterization**: `winding_number` (1D
  two-band), `chern_number` (Fukui–Hatsugai–Suzuki), `TKNN`
  (`σ_xy = C`), `ChernFromBerryCurvature` (`C = (1/2π)∫Ω`),
  `BulkBoundary` (`n = |ν|`, edge-mode count from the bulk invariant)
- **Heat-capacity relations**: `SpecificHeatFDT`, `SpecificHeatFromEntropy`
  (`c = T ∂s/∂T`), `HeatCapacityDifference` (Mayer `c_p − c_v = T v α²/κ_T`)
- **Thermodynamic identities**: the four `Maxwell*` relations (`Helmholtz`
  / `Gibbs` / `Internal` / `Enthalpy`), `ClausiusClapeyron`
  (`dp/dT = L/(T ΔV)`), `GibbsDuhem` (`S dT − V dp + N dμ = 0`)
- **Ensembles & thermal pure quantum states**: `MicrocanonicalTemperature`
  (`β = ∂S/∂E`, the microcanonical–canonical bridge), `CanonicalTPQ`
  (`Z = D·⟨ψ₀|e^{−βH}|ψ₀⟩`, Sugiura–Shimizu) — the identities an ED / MPS /
  TPQ calculation checks its measured quantities against
- **Structure-factor sum rules & dynamical FDT**: `StaticFromDynamicalStructureFactor`
  (`S(q) = ∫S(q,ω)dω/2π`, Van Hove), `DynamicalFDT`
  (`S(q,ω) = χ''/[π(1−e^{−βω})]`, reproduces detailed balance),
  `StructureFactorSusceptibility` (`χ = β S(q→0)`)
- **Entanglement entropy**: `RenyiTwoPurity` (`S_2 = −ln Tr ρ²`),
  `CFTEntanglementSlope` (`dS/d ln ℓ = c/3`, reads off the central charge;
  Calabrese–Cardy), `page_average_entropy` (Page's random-state formula)
- **CFT finite-size**: `CasimirCentralCharge` (`e₀(L)=e_∞−πcv/6L²`, reads off c), `FiniteSizeGap` (`E_x−E₀=2πvx/L`, reads off scaling dimensions) — extract universal data straight from finite-size MPS/ED spectra
- **Nonlinear-tensor symmetry**: `intrinsic_permutation_symmetric`,
- **Gaps & correlation length**: `CorrelationLengthGap` (`ξ = v/Δ`), `DynamicalScaling` (`d lnΔ/d lnξ = −z`, the dynamical critical exponent)
  `canonical_component`, `permutation_equivalent` — `χ⁽ⁿ⁾`'s field
  indices (with their frequencies) are interchangeable, so
  `χ⁽²⁾_{x;yz} = χ⁽²⁾_{x;zy}`
- **Dynamical / spectral identities**: `Dyson` (`G⁻¹ = G₀⁻¹ − Σ`),
  `SpectralFromGreens` (`A = −Im Gᴿ/π`), `SpectralSumRule` (`∫A dω = 1`),
  `DetailedBalance` (`S(q,−ω) = e^{−βω} S(q,ω)`), `NMRExponent`
  (`θ_NMR = 2Δ_op − 1`)

**Tensor structure — internal degrees of freedom.** Quantities that are
tensors are not blurred into scalars: they carry their indices as type
parameters and declare `tensor_rank` / `index_spaces` / `indices`.
`Susceptibility(:x, :y)` is the off-diagonal `χ_xy` the fused
`SusceptibilityZZ`-style names could not express, and the design is
**order-extensible to nonlinear response**: `Susceptibility(:x, :y, :z)`
is the second-order `χ⁽²⁾_{x;yz} = ∂²M_x/∂h_y∂h_z` (`response_order == 2`,
`tensor_rank == 3`), with the genealogy extending recursively `χ⁽ⁿ⁾ ⟵
χ⁽ⁿ⁻¹⁾ ⟵ … ⟵ M ⟵ F`. Index spaces: `SpinAxis` / `SpatialDirection`
(`Conductivity σ_μν`, also nonlinear) / `OrbitalIndex` (propagators
`G_ab`). `Dyson` is written with `inv`, so the one identity holds for
scalar single-band and matrix orbital-space propagators alike.

**Multi-time.** `frequency_arguments` is a quantity's multi-time
dimensionality — the number of independent frequency (⇔ time) variables.
The static `Susceptibility{I}` is the zero-frequency limit (`0`); the
**dynamical** `DynamicalSusceptibility(:x, :y, :z)` is
`χ⁽²⁾(ω₁, ω₂)` — an `n`-th order nonlinear response applies the field at
`n` distinct times, so `frequency_arguments == n` (2D coherent
spectroscopy). Its microscopic origin is the **Kubo formula** — the
retarded (`n`-fold nested-commutator) response function of the
correlation, `DynamicalSusceptibility ⟵ :kubo ⟵ DynamicalCorrelation`
in the dynamical graph. References for the response theory
([`Kubo1957`], [`WanArmitage2019`]) and the scaling / FDT relations
([`Rushbrooke1963`]…, [`CallenWelton1951`]) are in `docs/references.bib`,
each DOI CI-verified against Crossref.

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
- **Response genealogy**: `derivative_edge` encodes the derivative tree
  rooted at the free energy — `M = −∂F/∂h`, `χ = ∂M/∂h`, `C = ∂U/∂T`,
  `U ⟵ F`. `differentiation_chain(Susceptibility) ==
  [Susceptibility, Magnetization, FreeEnergy]`,
  `derivative_order(Susceptibility, MagneticField) == 2`
  (`χ = ∂²F/∂h²`). Exact formulas live in the paired relations
  (`MagnetizationResponse`, `SusceptibilityResponse`, `GibbsHelmholtz`);
  the *definitional* `χ = ∂M/∂h` and the *statistical* `χ = β·Var(M)`
  (`SusceptibilityFDT`) are the same response two ways.

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
