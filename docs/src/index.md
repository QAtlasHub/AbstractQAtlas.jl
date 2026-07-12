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

**Equalities and inequalities.** Most relations are equalities
(`check` ≡ `abs(residual) ≤ atol`); bound-type constraints are declared
with [`@inequality`](@ref) as [`AbstractInequality`](@ref), whose residual
is the `≥ 0` **slack** — `check` tests that direction, [`slack`](@ref)
reports the margin, and `solve` returns the **saturation** (tight-bound)
value.  The quantum-information entropy inequalities are the first users:
`EntropyNonNegativity`, `MaxEntropyBound` (`S ≤ ln d`), `Subadditivity`,
`ArakiLieb`, `StrongSubadditivity` (Lieb–Ruskai), `RenyiMonotonicity`.

```julia
check(StrongSubadditivity(); S_AB, S_BC, S_ABC, S_B)   # S_AB + S_BC ≥ S_ABC + S_B ?
solve(Subadditivity(), Val(:S_AB); S_A, S_B)           # the tight bound S_A + S_B
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

## Tensor structure — internal degrees of freedom

Quantities that are tensors carry their indices as type parameters and
declare [`tensor_rank`](@ref) / [`index_spaces`](@ref) / [`indices`](@ref),
so they are not silently scalarized. `Susceptibility(:x, :y)` is the
off-diagonal `χ_xy`; the design is **order-extensible to nonlinear
response** — `Susceptibility(:x, :y, :z)` is `χ⁽²⁾_{x;yz} = ∂²M_x/∂h_y∂h_z`
([`response_order`](@ref) `== 2`, [`tensor_rank`](@ref) `== 3`) and the
[`differentiation_chain`](@ref) extends recursively `χ⁽ⁿ⁾ ⟵ χ⁽ⁿ⁻¹⁾ ⟵ … ⟵ M
⟵ F`. Index spaces: [`SpinAxis`](@ref), [`SpatialDirection`](@ref)
([`Conductivity`](@ref), also nonlinear), [`OrbitalIndex`](@ref)
(propagators). [`Dyson`](@ref) is written with `inv`, so the identity
holds for scalar single-band and matrix orbital-space propagators alike.

[`frequency_arguments`](@ref) is a quantity's **multi-time
dimensionality** — the count of independent frequency (⇔ time)
variables. The static [`Susceptibility`](@ref) is the zero-frequency
limit (`0`); the dynamical [`DynamicalSusceptibility`](@ref)`(:x, :y, :z)`
is `χ⁽²⁾(ω₁, ω₂)` with `frequency_arguments == 2` — an `n`-th order
nonlinear response is intrinsically multi-time (2D coherent
spectroscopy). Its microscopic origin is the **Kubo formula**
([`spectral_origin`](@ref)`(DynamicalSusceptibility(:x,:y,:z)) ==
(DynamicalCorrelation{(:x,:y,:z)}, :kubo)`): the retarded `n`-fold
nested-commutator response function of the **same-order** correlation —
an `n`-th order response is an `n`-time (`(n+1)`-point) correlation, so
the Kubo edge preserves the frequency count on both sides.
Response-theory and scaling/FDT references (Kubo 1957,
Wan–Armitage 2019, Rushbrooke/Widom/Fisher/Josephson, Callen–Welton
1951) live in `docs/references.bib`, DOI-verified in CI.

## Nonlinear-tensor symmetry & accumulated relations

The nonlinear susceptibility is *essentially* a higher-order tensor, and
carries **intrinsic permutation symmetry** — `χ⁽ⁿ⁾`'s field indices (with
their frequencies) are interchangeable, so
`Susceptibility(:x, :y, :z) == Susceptibility(:x, :z, :y)` under
[`permutation_equivalent`](@ref) ([`canonical_component`](@ref) sorts the
field indices; the response index is fixed).

Known inter-quantity relationships accumulate as first-class relations:
`ChernFromBerryCurvature` (`C = (1/2π)∫Ω`) with [`TKNN`](@ref)
(`σ_xy = C`) and `BulkBoundary` (`n = |ν|`) for the topological side;
`SpecificHeatFromEntropy` (`c = T ∂s/∂T`) and `HeatCapacityDifference`
(Mayer's `c_p − c_v = T v α²/κ_T`) for heat capacity; `MicrocanonicalTemperature`
(`β = ∂S/∂E`, the microcanonical–canonical bridge) and `CanonicalTPQ`
(`Z = D·⟨ψ₀|e^{−βH}|ψ₀⟩`, Sugiura–Shimizu) for statistical ensembles and
thermal-pure-quantum estimators.  The **transport** family
([`Conductivity`](@ref) with its AC [`DynamicalConductivity`](@ref),
[`ThermalConductivity`](@ref), [`Thermopower`](@ref),
[`PeltierCoefficient`](@ref), [`DrudeWeight`](@ref), and the
[`ElectricCurrent`](@ref) / [`HeatCurrent`](@ref)) carries
`WiedemannFranz` (`κ = L₀σT`), the Mott `MottFormula`
(`S = −(π²/3)T d ln σ/dε`), the Kelvin `KelvinRelation` (`Π = TS`),
`OnsagerReciprocity` (`L_{μν} = L_{νμ}`), and the optical `OpticalSumRule`
(`∫Re σ dω = πD + W_reg`).  References are DOI-verified in
`docs/references.bib`.

**Scope note.** The Berry curvature is the *imaginary* part of the
quantum geometric tensor; the real part (the quantum metric) and the
mixed-state / Uhlmann generalizations are deliberately **out of scope** —
this package stays at the model-independent textbook level.

## Structure — definitional correspondences

The `structure/` layer holds the generic facts that are true *by
definition*, from which the forms above are derived rather than
restated:

- **Critical correspondence** — [`critical_scaling`](@ref) maps each
  observable to the exponent governing its singularity
  (`Susceptibility ↦ γ`, `SpontaneousMagnetization ↦ β`, …). The
  singular and finite-size forms follow: [`singular_form`](@ref),
  [`fss_size_exponent`](@ref) (so `χ_max ∼ L^{γ/ν}` is *derived*, not a
  hand-passed ratio), [`fss_peak`](@ref), [`collapse_coordinates`](@ref).
  Field-driven and distance-driven laws: [`critical_isotherm`](@ref)
  (δ), [`correlation_decay`](@ref) (η).
- **Transition classification** — [`FirstOrder`](@ref) /
  [`ContinuousTransition`](@ref) / [`KosterlitzThouless`](@ref), each
  carrying [`ehrenfest_order`](@ref), [`has_order_parameter`](@ref),
  [`has_latent_heat`](@ref), [`has_critical_exponents`](@ref).
- **Response genealogy** — [`derivative_edge`](@ref) encodes the
  derivative tree rooted at the free energy (`M = −∂F/∂h`, `χ = ∂M/∂h`,
  `C = ∂U/∂T`). [`differentiation_chain`](@ref) traces any response back
  to [`FreeEnergy`](@ref); [`derivative_order`](@ref) counts field
  derivatives (`χ = ∂²F/∂h²` ⇒ order 2); [`potential_root`](@ref),
  [`is_response`](@ref), [`conjugate_field`](@ref). The exact formulas
  are the paired relations [`MagnetizationResponse`](@ref),
  [`SusceptibilityResponse`](@ref), [`GibbsHelmholtz`](@ref) — the
  definitional `χ = ∂M/∂h` and the statistical `χ = β·Var(M)`
  ([`SusceptibilityFDT`](@ref)) being the same response two ways.
- **Dynamical / spectral graph** — [`spectral_origin`](@ref) /
  [`spectral_chain`](@ref) trace the frequency-resolved quantities back
  to their sources (`DensityOfStates ⟵ SpectralFunction ⟵
  RetardedGreensFunction ⟵ SelfEnergy`), and [`origin_relation`](@ref)
  ties each single-`(q,ω)`-point edge to its exact identity: [`Dyson`](@ref)
  (`G⁻¹ = G₀⁻¹ − Σ`), [`SpectralFromGreens`](@ref) (`A = −Im Gᴿ/π`),
  with [`SpectralSumRule`](@ref), [`DetailedBalance`](@ref)
  (`S(q,−ω) = e^{−βω} S(q,ω)`) and [`NMRExponent`](@ref)
  (`θ_NMR = 2Δ_op − 1`). Transform / BZ-sum / limit edges have no
  pointwise form — their *evaluation* is the functional sibling's job.

## API reference

```@autodocs
Modules = [AbstractQAtlas]
```
