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

## The three-verb interface

Every relation implements a uniform contract:

- `residual(rel; vars...)` — signed violation; `0` ⇔ satisfied,
- `check(rel; atol=0, vars...)` — `|residual| ≤ atol`,
- `solve(rel, Val(:x); vars...)` — the value of `x` implied by the rest.

with an **exact-arithmetic contract**: `Rational` in ⇒ `Rational` out,
so exactly-known values satisfy their relations exactly, not merely to
floating-point tolerance.

```julia
using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

# 2D Ising exponents are exact rationals — the relations hold EXACTLY:
residual(Rushbrooke(); α=0//1, β=1//8, γ=7//4)   # 0//1
solve(Widom(), Val(:γ); β=1//8, δ=15//1)          # 7//4
check(Fisher(); γ=7//4, ν=1//1, η=1//4)           # true

# gate a whole exponent table at once:
exponents_consistent((α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4); d=2)  # true

# fluctuation–dissipation: solve() doubles as the estimator
# c_v = β² Var(E) / N
# solve(SpecificHeatFDT(), Val(:C); var_E=var_E, β=β, N=N)
```

## What is covered (v0.1)

- **Scaling laws** — [`Rushbrooke`](@ref), [`Widom`](@ref),
  [`Fisher`](@ref), [`Josephson`](@ref), gated by
  [`exponents_consistent`](@ref).
- **Fundamental equations** — [`FreeEnergyFromZ`](@ref)
  (``F = -β^{-1}\ln Z``), [`FreeEnergyLegendre`](@ref) (``F = U - TS``),
  [`EntropyResponse`](@ref) (``S = -∂F/∂T``), [`GibbsHelmholtz`](@ref)
  (``U = ∂(βF)/∂β``).
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
