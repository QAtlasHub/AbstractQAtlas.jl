# Design: the core function surface

Status: **establishing** — this document fixes the canonical *core-function*
design of AbstractQAtlas (the verify-engine's public verb layer) so it stops
drifting as relations accumulate. It is the source of truth for the five design
pillars below. Pillars 1–2 mostly *codify* behaviour that already exists (the
deliverable is this document plus small hardening); pillars 3–5 name concrete
gaps that follow-on feature PRs close. It is the sibling of
[`type-keyed-interface.md`](type-keyed-interface.md), which this generalises
from the relation object to the whole verb surface.

## 0. What "core" is — and what is not

The package exports ~90 names. They fall into two layers, and only the first is
*core*:

- **Core** — the verify-engine's fundamental verbs (this document). Stable
  public API; changes are semver-significant.
- **Domain content** — the physics accessors (`critical_exponent`,
  `wick_pfaffian`, `spectral_origin`, `fss_peak`, `keldysh_distribution`, …).
  These are a relation family's *own vocabulary*, added and removed as the
  physics accumulates; they are not the engine.

The core is exactly three surfaces:

```
A. the RELATION        @relation/@inequality  ·  residual/check/solve/slack  ·  introspection
B. ADOPTION            bag  →  applicable_relations / relation_report / check_all
C. the SEAMS           fetch · report · principal_value_hilbert/spectral_moment · thermal_derivative · derive
                       ("interface here, values / numerics / derivations there")
```

Everything below is one of these three. The five pillars are the design
principles that make them *certain*, *extensible*, and *trustworthy*.

---

## Pillar 1 — Type-based certainty: the type IS the identity

**Principle.** A relation's variables are identified by their **quantity / field
type**, never by a formula-letter symbol. The type is collision-free and
drift-free; the symbol (`G`, `A`, `S`) is a *private local binding* for the
math in the kernel, not an external contract.

**Why (the correctness hazard, not ergonomics).** A symbol is chosen
independently of the type, so it *drifts* (`RetardedGreensFunction` is `G` in
`Dyson` but `GR` in the Keldysh relations) and *collides* (`:S` is four distinct
quantities — thermal entropy, thermopower, dynamical structure factor, von
Neumann entropy). On the symbol path `derive(:Π; S = entropy, T)` once returned
`6.0`, silently reusing entropy as the Seebeck coefficient — the verify-engine
producing a **false success inside itself**. Types make that impossible.

**Established contract.**

- **Declaration.** `@relation :domain Name(v::Type, …) = expr`. The `::Type`
  annotation is the **bag key (metadata)**, not a value constraint — a
  matrix-valued `Dyson` still type-checks. An *untyped* slot is a **supplied
  scalar** (evaluation coordinates `ω`/`q`, supplied integrals `pv_imag`/`∫A`,
  distribution factors `h`/`F`): identity-bearing vs supplied-slot is the split
  from [`type-keyed-interface.md`](type-keyed-interface.md) §R3.
- **The canonical carrier is the `Bag`.** `bag(Type => value, …)` is the
  semantic ground truth; the engine-level verbs (`relation_report`,
  `check_all`, `applicable_relations`, `derive`, `derivable`) are **bag-first**
  and key on `VariableKey(type, support)`.
- **Symbol-kwargs is an ergonomic shortcut, not the ground truth.**
  `residual(rel; G = …, Σ = …)` stays for convenient single-relation calls, but
  it is explicitly *not* how the engine reasons about identity — it is the form
  that drifts and collides, tolerated only where a human names one relation's
  variables directly.
- **Load-time guards make a mis-declaration fail at load, not at runtime.** The
  macro's `_validate_relation` requires identity types to be **concrete** (a
  bare `Susceptibility` family cannot `===`-match a concrete bag key → silent
  dead relation) and **distinct** (two slots of the same type bind the same
  value → an always-pass relation, the engine's own worst failure), and allows
  **≤ 1** `UnionAll` family slot (family auto-discovery, §8a of the sibling doc).
- **Certainty reaches the graph.** `derive` / `typed_derivation_graph` operate
  over `VariableKey` nodes, so the `:S`-collision derivation bug is
  *structurally* absent, not merely untriggered.

**Gap / hardening.** The undecorated + tensor + transport domains are migrated
(#79–#85); the decorated spectral `±ω`/Re–Im and the entanglement Region epic
still carry symbol keys by design (they need the support/decoration machinery).
No new work is *required* for the contract to hold — it is the established
default; new relations are type-keyed unless they have no named subject
(generic laws: Maxwell, Onsager, Ehrenfest stay symbol-keyed on purpose).

---

## Pillar 2 — The `fetch` seam: one non-exported retrieval verb, everything routes through it

**Principle.** A stored/computed *value* enters the engine through exactly one
canonical verb, `fetch(model, quantity, bc)`. There is no second retrieval
path.

**Established contract.**

- **Generic function + informative fallback.** `fetch` is owned by
  AbstractQAtlas as a generic function whose top-level method **errors
  informatively** ("no fetch method for model=…, quantity=…, bc=…; the
  implementing package must define …"). Each implementing package (QAtlas)
  registers **one method per `(model, quantity, bc)` triple**.
- **Deliberately NOT exported.** `fetch` clashes with `Base.fetch`, so it is
  called qualified (`AbstractQAtlas.fetch`) or imported explicitly
  (`using AbstractQAtlas: fetch`). This is the **canonical convention**, not a
  wart — QAtlas does exactly the same. New consumers follow it; we do **not**
  add a second non-clashing alias (one true name).
- **The seam pattern is the ecosystem's single extension mechanism.** *Generic
  verb + erroring fallback in AbstractQAtlas; methods at the leaves.* Every seam
  instantiates it, so the package is the `AbstractFFTs` of the atlas family:

  | seam | verb(s) | who implements | carries |
  |---|---|---|---|
  | reference values | `fetch(model, quantity, bc)` | QAtlas | stored/oracle values |
  | reported values | `report(model, quantity, bc; …) → Card` | any reporter | a computed value + provenance |
  | functional numerics | `principal_value_hilbert` / `spectral_moment(::AbstractResponse, …)` | functional sibling | transforms / quadratures |
  | autodiff | `thermal_derivative(quantity, potential, x)` | ForwardDiff ext | genealogy derivatives (Pillar 4) |

- **Cycle-free.** Consumers depend on AbstractQAtlas *only* (never the reverse);
  the values/numerics live at the leaves. `report` is the reporter-facing
  sibling of `fetch` — reporters push DATA (a `Card`), they are not depended on.
- **"Route through `fetch`."** Higher verbs that need a value — a `derive`
  cross-check, a round-trip against a reported `Card` — obtain it via `fetch`,
  never through a parallel lookup. One place a value enters ⇒ one place to
  cache (Pillar 5) and one place to trust.

**Gap / hardening.** The contract holds today; the remaining work is to *state*
it (this document) and to keep new seams (e.g. any future graph/data-export
verb) inside the pattern rather than inventing a parallel mechanism.

---

## Pillar 3 — Citation discipline: every relation carries a verified reference

**Principle.** A verify-engine that cannot cite its own laws cannot be trusted.
Every physics relation carries a reference, and every reference is real (ties
the `never-fabricate-citations` rule).

**Established contract.**

- **Every `@relation`/`@inequality` docstring carries a reference**, in one of
  two forms:
  - **`[key](@cite)`** to a `docs/references.bib` entry, for any result with a
    DOI or arXiv id. The entry is **doiget-verified** before it is added.
  - **Honest inline** ("Langreth (1976)"; "Haug & Jauho, *Quantum Kinetics …*";
    "Zamolodchikov, JETP Lett. 43, 730 (1986)"), for pre-DOI / textbook / book
    results with no canonical resolvable id — **no fabricated bib entry**.
- **Never invent a real-looking DOI.** The `citations` CI gate (the doiget
  `verify` action) checks every bib DOI/arXiv resolves on Crossref/arXiv, so a
  fabricated or mistyped id fails CI. A real-but-unresolvable id is excepted by
  hand in `docs/references.allow`.
- **Enforcement (new).** A test asserts that **every registered relation's
  docstring contains a reference token** — a `[…](@cite)` *or* an inline
  `(… YYYY)` year pattern — so a newly-added uncited relation fails CI, exactly
  as the reflection invariants guard tensor traits.

**Gap / roadmap.** `thermodynamic.jl` (16 relations) and `fundamental.jl` (7)
carry **no `@cite`** — the textbook thermodynamics (Maxwell, Clausius–Clapeyron,
Gibbs–Duhem) and the `Z`–`F`–`U`–`S` web. Audit each: attach an honest inline
textbook citation where no single canonical paper exists, and a doiget-verified
bib entry where one does; then land the enforcement test. This is the most
self-contained pillar and the recommended first execution PR.

---

## Pillar 4 — The AD extension seam: derivatives from the response genealogy

**Principle.** The response genealogy **declares** that a quantity is a
derivative of a potential (`derivative_edge`: `M = −∂F/∂h`, `χ = ∂M/∂h`,
`C = ∂U/∂T`, …); automatic differentiation **computes the value**. Interface in
core, backend in an extension.

**Established contract.**

- **Generic function + optional-backend extension.** `thermal_derivative(quantity,
  potential, x)` is a stub with an erroring fallback in `src/autodiff.jl`; the
  methods live in the **ForwardDiff package extension**
  (`ext/AbstractQAtlasForwardDiffExt.jl`), so AD is an *optional* dependency.
- **Backend-agnostic by construction.** `residual`/`check`/`solve` and the
  scaling/response *forms* are pure arithmetic and already differentiate through
  *any* AD backend with no extension — the ext adds only the genealogy
  AD-*evaluation* (potential **function** → derived-quantity **value**). A
  reverse-mode backend is therefore a *second* extension on the same generic
  function, not a redesign.
- **Genealogy-structured, order-exact.** The derivative order and field are read
  from `response_order`/`derivative_edge`, so
  `thermal_derivative(Susceptibility(α, β₁…βₙ), F, h⃗, components)` is the exact
  mixed partial `−∂ⁿ⁺¹F/∂h_α∂h_{β₁}…`; a single-field `F(h)` fixes only the
  **diagonal** and *errors* on an off-diagonal request rather than silently
  returning the diagonal.

**Current coverage.** `Magnetization`, `Susceptibility` (diagonal + multi-field
mixed partial), `ThermalEntropy`, `SpecificHeat`, `Energy` (Gibbs–Helmholtz).

**Gap / roadmap.** Extend along the remaining genealogy edges — the dynamical
response (frequency derivatives of `DynamicalSusceptibility`), transport
coefficients as field-derivatives — driving each from `derivative_edge` so the
map is not hand-maintained. Optionally add a reverse-mode ext (Zygote/Enzyme)
for many-parameter potentials: same verb, second `[extensions]` entry.

---

## Pillar 5 — Efficient quantity retrieval: cache the pure lookups

**Principle.** A value keyed by `(model, quantity, bc)` is a **pure** function of
that triple (the reference-value contract already implies determinism), and the
oracles behind it (ED / DMRG / TPQ) are expensive. Computing once and caching is
therefore safe and often necessary.

**Established contract.**

- **Retrieval efficiency is a `fetch`-layer concern (Pillar 2).** Because every
  value enters through `fetch`, the cache keys on the fetch triple
  `(model, quantity, bc)` (plus the size/`N` kwarg where a method reads it) and
  memoizes the pure result — one seam, one cache.
- **Opt-in and purity-explicit.** Caching assumes a `fetch` method is a
  deterministic pure function of its triple; a method with side effects (or a
  live-recomputed value) opts out. The cache is invalidated only by an explicit
  reset — no time-based expiry, no hidden staleness.
- **stdlib-only.** An `IdDict`/`Dict` in the dispatch layer (or a `fetch_cached`
  wrapper), matching the existing `_DERIV_STEPS`/`_TYPED_DERIV_STEPS`
  built-once-and-cached discipline in `derivation.jl` — no `Memoize.jl`
  dependency.
- **`derive` inherits it.** A *derived* quantity's route and value memoize the
  same way, so a repeated `derive(target, bag)` does not re-walk the graph.

**Gap / roadmap.** No cache exists today (each `fetch` re-dispatches; each
`derive` re-solves). The concrete PR adds the memoization layer *after* the
`fetch` seam contract (Pillar 2) is codified, so the cache sits at the one
canonical entry point.

---

## Roadmap — execution order

1. **This document** — establish the five contracts. (Pillars 1–2 are largely
   codification; their deliverable is mostly here.)
2. **Citation audit + enforcement test (Pillar 3)** — self-contained, immediate;
   fill `thermodynamic.jl`/`fundamental.jl`, then land the "every relation
   cites" guard.
3. **AD extension (Pillar 4)** — extend `thermal_derivative` along the genealogy;
   optional reverse-mode ext.
4. **`fetch` cache (Pillar 5)** — the memoization layer at the `fetch` seam,
   after Pillar 2 is fixed.

Each is a focused PR; pillars 1–2 additionally get small hardening PRs only
where the code and this contract disagree.
