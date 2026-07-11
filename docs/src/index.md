# AbstractQAtlas.jl

*The model-independent layer of the QAtlas ecosystem — abstract quantity
vocabulary + generic physics relations as first-class, tested objects.*

In the spirit of `AbstractFFTs`: concrete atlases
([QAtlas.jl](https://github.com/sotashimozono/QAtlas.jl)) **implement**
this package, never the reverse.

## Division of responsibility

| lives here (AbstractQAtlas) | lives in the implementing atlas |
|---|---|
| type vocabulary: `AbstractQAtlasModel`, `AbstractQuantity`, `BoundaryCondition`, the generic `fetch` verb | concrete models and registered `fetch` methods |
| generic relations: scaling laws, fluctuation–dissipation, Wick's theorem, topological invariants, FSS forms | reference **values** (critical temperatures, exact magnetizations, exponent tables) |

A relation is an *identity among observables or exponents* — a
statement true independently of any model.  Expressing each one once,
as a tested object, means downstream packages stop re-deriving them ad
hoc in comments and per-model tests.

## Declare once, derive everything

A relation is written exactly once, with [`@relation`](@ref):

```julia
@relation :scaling Rushbrooke(α, β, γ) = α + 2β + γ - 2
```

One declaration yields the struct, the residual kernel, the
[`variables`](@ref)/[`domain`](@ref) introspection traits, registry
membership, and — with no hand-written rearrangements — [`solve`](@ref)
for every variable the expression is affine in (a non-affine variable is
*refused*, never silently mis-solved).

The uniform verbs:

- [`residual`](@ref)`(rel; vars...)` — signed violation; `0` ⇔ satisfied,
- [`check`](@ref)`(rel; atol=0, vars...)` — `|residual| ≤ atol`,
- [`solve`](@ref)`(rel, Val(:x); vars...)` — the value of `x` implied by the rest,

with an **exact-arithmetic contract**: `Rational` in ⇒ `Rational` out,
so exactly-known values satisfy their relations exactly, not merely to
floating-point tolerance.  Relations taking an inverse temperature
accept `β` or `T` at every verb; normalization happens once, in the verb
layer.

```julia
using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

residual(Rushbrooke(); α=0//1, β=1//8, γ=7//4)   # 0//1 — exact
solve(Widom(), Val(:δ); β=1//8, γ=7//4)           # 15//1 — derived, not hand-coded
check(Fisher(); γ=7//4, ν=1//1, η=1//4)           # true
```

## Adopting from another package: one call

A consumer never hand-lists relations — [`applicable_relations`](@ref)
selects by variable names, [`relation_report`](@ref) evaluates and
reports, [`check_all`](@ref) gates (an empty match is `false`, never a
silent green):

```julia
check_all((α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4, d=2))   # exponent table gate
relation_report((C=c, var_E=v, T=T, N=N); atol=tol)                  # thermodynamics sweep
```

Pass `domain=` when a data set mixes families (physics overloads names:
the exponent `β` vs the inverse temperature `β`).  Downstream packages
declare their own relations with the same [`@relation`](@ref) macro.

## What is covered (v0.1)

- **Scaling laws** — [`Rushbrooke`](@ref), [`Widom`](@ref),
  [`Fisher`](@ref), [`Josephson`](@ref), gated by
  [`exponents_consistent`](@ref).
- **Fundamental equations** — [`FreeEnergyFromZ`](@ref)
  (``F = -β^{-1}\ln Z``), [`FreeEnergyLegendre`](@ref) (``F = U - TS``),
  [`EntropyResponse`](@ref) (``S = -∂F/∂T``), [`GibbsHelmholtz`](@ref)
  (``U = ∂(βF)/∂β``).
- **Distributions & statistics** — [`MicroCanonical`](@ref) /
  [`Canonical`](@ref) / [`GrandCanonical`](@ref) / [`Squeezed`](@ref)
  with [`ensemble_weight`](@ref); [`Fermionic`](@ref)/[`Bosonic`](@ref)
  [`occupation`](@ref) functions; the [`ThermalAverage`](@ref) marker
  (``⟨Q⟩_D`` at the type level).
- **Thermodynamic identities** — [`SpecificHeatFDT`](@ref),
  [`SusceptibilityFDT`](@ref), [`LinearResponseFDT`](@ref)
  (conventions: ``c_v = β^2\,\mathrm{Var}(E)/N``,
  ``χ = β\,\mathrm{Var}(M)/N``).
- **Wick's theorem** — [`wick_contraction`](@ref),
  [`wick_density_correlation`](@ref) (number-conserving Gaussian
  fermion states; BdG/Pfaffian generalization tracked).
- **Topological invariants** — [`winding_number`](@ref) (1D two-band),
  [`chern_number`](@ref) (Fukui–Hatsugai–Suzuki), [`TKNN`](@ref).
- **FSS forms** — [`collapse_coordinates`](@ref),
  [`fss_peak_scaling`](@ref) and friends: the vocabulary of a
  finite-size-scaling analysis.

## API reference

```@autodocs
Modules = [AbstractQAtlas]
```
