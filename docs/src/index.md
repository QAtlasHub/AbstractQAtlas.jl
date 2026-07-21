# AbstractQAtlas.jl

*The model-independent layer of the QAtlas ecosystem — abstract quantity
vocabulary + generic physics relations as first-class, tested objects.*

In the spirit of `AbstractFFTs`: concrete atlases
([QAtlas.jl](https://github.com/QAtlasHub/QAtlas.jl)) **implement**
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
- **The exponents' RG origin** — the four scaling laws
  ([`Rushbrooke`](@ref), [`Widom`](@ref), [`Fisher`](@ref),
  [`Josephson`](@ref)) are not four independent axioms: they all follow
  from the homogeneity of the singular free energy
  `f_s(t,h) = b^{-d} f_s(b^{y_t}t, b^{y_h}h)`, i.e. from just two RG
  eigenvalues plus `d`. [`ScalingDimensions`](@ref)`(y_t, y_h, d)` carries
  them and [`critical_exponents`](@ref) *derives* the whole set
  `(α,β,γ,δ,ν,η)` — for which every scaling relation has residual `≡ 0` by
  construction ([`exponents_consistent`](@ref) is `true` for any input).
  The inverse [`scaling_dimensions`](@ref)`(; ν, η, d)` recovers the
  eigenvalues, so `δ`, `β`, `γ`, `α` are all reconstructible from `(ν, η,
  d)` alone. Values are inputs; the map is the universal content:

  ```julia
  critical_exponents(ScalingDimensions(1//1, 15//8, 2))
  # (α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4)  ← 2D Ising, exact
  ```
- **Transition classification** — [`FirstOrder`](@ref) /
  [`ContinuousTransition`](@ref) / [`KosterlitzThouless`](@ref), each
  carrying [`ehrenfest_order`](@ref), [`has_order_parameter`](@ref),
  [`has_latent_heat`](@ref), [`has_critical_exponents`](@ref).
- **Response genealogy** — [`derivative_edge`](@ref) encodes the
  derivative tree rooted at a thermodynamic potential (`M = −∂F/∂h`,
  `χ = ∂M/∂h`, `C = ∂U/∂T` at the free-energy root; `N = −∂Ω/∂μ` at the
  grand-potential root). [`differentiation_chain`](@ref) traces any response
  back to its [`potential_root`](@ref) ([`FreeEnergy`](@ref) or
  [`GrandPotential`](@ref)); [`derivative_order`](@ref) counts field
  derivatives (`χ = ∂²F/∂h²` ⇒ order 2); [`is_response`](@ref),
  [`conjugate_field`](@ref). The exact formulas
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
- **Keldysh RAK structure & fluctuation–dissipation** — the real-time
  sibling of the Matsubara Dyson layer. The retarded/advanced/Keldysh
  components obey the state-independent identities `G^K = G^> + G^<` and
  `G^R − G^A = G^> − G^<` ([`KeldyshComponent`](@ref),
  [`KeldyshCausality`](@ref)) plus the adjoint tie `G^A = (G^R)†`
  ([`AdvancedRetardedConjugate`](@ref), which also puts `G^A` on the
  spectral graph via the `:adjoint` edge). In equilibrium the Keldysh
  component is *not* independent: the **fluctuation–dissipation theorem**
  [`KeldyshFDT`](@ref) locks it to the spectral part,
  `G^K = h(ω)(G^R − G^A)`, with the distribution function
  [`keldysh_distribution`](@ref) `h = coth(βω/2)` (bosons) `= tanh(βω/2)`
  (fermions) `= 1 ∓ 2n`. That `h` is itself a *consequence* of the
  [`KMSGreaterLesser`](@ref) detailed-balance relation
  `G^<(ω) = ζe^{−βω}G^>(ω)` — FDT is derived, not assumed — and
  [`SpectralFromKeldysh`](@ref) `A = i(G^R−G^A)/2π` reduces to `−Im G^R/π`,
  bridging back to the normalized spectral function.

## One queryable graph of physics

The genealogy, the spectral graph, the Fourier pairs and the
relation ↔ quantity links are all the same shape — a *typed edge between
quantity kinds* — so they fold into **one** queryable graph
([`quantity_graph`](@ref)), mirroring the vocabulary of QAtlas's *model*
graph (`relations(model)`): the two atlases share one graph language,
models ⊕ quantities.

- [`related_quantities`](@ref)`(q)` — the neighborhood of a quantity: every
  edge it participates in, tagged by kind (`:derivative`, `:spectral`,
  `:fourier`, `:law`) as a [`QuantityEdge`](@ref). For `Susceptibility` this
  surfaces both `Magnetization` (as `∂/∂h` **and** via the FDT) and
  `StaticStructureFactor` (via the structure-factor sum rule) at once.
- [`quantity_path`](@ref)`(a, b)` — the machine answer to *"how are `a` and
  `b` related?"*: a shortest path of typed edges. `SpecificHeat` and
  `Magnetization` connect through their shared [`FreeEnergy`](@ref) root.
- [`quantity_neighbors`](@ref)`(fam)` — incident edges in both directions
  (so a derivative *root* like `FreeEnergy` still surfaces the quantities
  that point at it), and [`quantity_graph_jsonl`](@ref) streams the whole
  network as stdlib-only JSONL for a graph view.

Nodes are quantity *families* (the index-erased `Susceptibility`, not
`Susceptibility{(:z,:z)}`) so the structural graph is finite; the concrete
index is still used to *resolve* an edge (`χ⁽²⁾ ⟶ χ⁽¹⁾ ⟶ M`) before the
endpoints collapse to their families.

Under the hood this is one instance of a generic parent, the
[`KnowledgeGraph`](@ref)`{N}` kernel — a bag of [`TypedEdge`](@ref)s over nodes
of type `N`, with the traversal/reachability/shortest-path/export
([`graph_reachable`](@ref), [`graph_shortest_path`](@ref),
[`graph_jsonl`](@ref)) written once.  `quantity_graph()` is a
`KnowledgeGraph{Type}`; the derivation graph below is a `KnowledgeGraph{Symbol}`
([`derivation_graph`](@ref)); and QAtlas's model graph becomes a third instance
once refactored onto this package — so a single graph-view renders
models ⊕ quantities ⊕ derivations.

## Deriving a quantity by route

Read each equality as a computation — [`solve`](@ref) turns a relation into
"given all-but-one variable, produce the last" — and the whole registry
becomes a **directed derivation graph** ([`derivation_steps`](@ref)): a node
per variable, a directed edge for every variable a relation can be solved for.
From a set of known quantities you can then ask what else is reachable and,
lazily, get one route run for you:

```julia
derivable(; Z = 2.0, β = 1.0)                # Set([:Z, :β, :f, …])  reachable
derive(:f; Z = 2.0, β = 1.0)                 # -0.6931…   (F = −β⁻¹ ln Z)
derive(:δ; α = 0//1, β = 1//8)               # 15//1  — two exact hops, Rushbrooke→Widom
```

Two safety guarantees, because a *chained* result is weaker than a directly
implemented one:

- **Equalities only.** Inequalities are excluded — their `solve` returns a
  *saturation bound*, not an equational value, so a bound can never
  masquerade as a derived quantity.
- **Auditable provenance.** [`derive`](@ref)`(…; debug=true)` returns a
  [`DerivationTrace`](@ref) — the value, the exact route (which relation
  produced each intermediate, from which inputs), and an `indirect` flag —
  so an indirectly derived number is inspected by its route, not trusted
  blindly. The route is discovered by calling the *real* `solve`, so a step
  whose relation is non-affine in its output is skipped, never faked; an
  unreachable target fails loudly.

```julia
derive(:δ; α = 0//1, β = 1//8, debug = true)
# DerivationTrace(:δ = 15//1  [indirect])
#   1. Rushbrooke: {α, β} → :γ
#   2. Widom: {β, γ} → :δ
```

## The scope line: definitional vs functional

Where does a dynamical quantity's *value* come from — this package or the
future ParaLA-based functional sibling? [`operation_scope`](@ref) draws
the line (issue #14):

- **`:definitional`** (here) — a **pointwise** identity relating quantity
  *values* at a single `(q, ω)`, or a *supplied* scalar (an integral, a
  derivative): `Dyson`, `SpectralFromGreens`, the sum rules,
  `KramersKronigReal`/`Imag`, every `@relation`. This package holds these
  as stdlib-only tested identities.
- **`:functional`** (sibling) — a **transform / sum / limit** that
  represents a quantity as a *function* and acts on it globally: the BZ
  average, the space-time Fourier transform, an `ω → 0` limit, the Kubo
  response. Only the structural edge lives here; the numerical evaluation
  is deferred.

The **grey zone** (cf. #6) resolves by the *supplied-value* convention: a
sum rule or Kramers–Kronig relation is `:definitional` — the relation
*checks a supplied number* here, while *computing* that number from the
function (the principal-value Hilbert transform, the spectral integral) is
`:functional`, the sibling's job. So the boundary is exactly
[`origin_relation`](@ref)'s split: definitional ⟺ a pointwise `@relation`
exists.

## API reference

```@autodocs
Modules = [
    AbstractQAtlas,
    AbstractQAtlas.StatisticalMechanics,
    AbstractQAtlas.Criticality,
    AbstractQAtlas.Correlations,
    AbstractQAtlas.Transport,
    AbstractQAtlas.QuantumInformation,
    AbstractQAtlas.QuantumFoundations,
    AbstractQAtlas.Topology,
]
```
