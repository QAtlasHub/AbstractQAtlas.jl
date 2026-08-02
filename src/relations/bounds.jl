# relations/bounds.jl — the universal bounds, stated against their bounding values.
#
# Every bound here has the same shape, and it is the shape `@bound`'s comparison
# form exists for: the BOUNDING side is a quantity (`CHSHBound`, `ChaosBound`, …)
# whose value a theory fixes and an atlas can fetch, while the BOUNDED side is an
# untyped slot — the measured correlator, the extracted Lyapunov exponent, the
# orthogonalization time someone timed.  No quantity in this vocabulary names
# those, and inventing one would be naming a measurement protocol, not a
# physical quantity.  Stating the bound anyway is the point: the relation is
# declared and discoverable now, and the untyped slot is filled by whoever has
# the measurement.  See docs/design/bounds.md §4.
#
# This also means none of these can be auto-materialized by a consumer that
# fetches both slots — by construction, not by omission.  A consumer supplies
# the bounded value itself.
#
# Domain tags: `:quantum` for the seven quantum-mechanical bounds, `:holographic`
# for Bekenstein, whose statement is gravitational.
#
# References (doiget-verified, docs/references.bib): [CHSH1969](@cite),
# [Tsirelson1980](@cite), [PopescuRohrlich1994](@cite), [Mermin1990](@cite),
# [MaldacenaShenkerStanford2016](@cite), [Bekenstein1981](@cite),
# [MargolusLevitin1998](@cite), [SekinoSusskind2008](@cite),
# [ShorPreskill2000](@cite), [BuzekHillery1996](@cite).

# ─── Bell-type correlator bounds ────────────────────────────────────────

"""
    CHSHInequality <: AbstractInequality

The CHSH inequality ([CHSH1969](@cite)): a measured CHSH correlator cannot
exceed what the theory admits,

`S ≤ S_max`

(slack `S_max − S`).  Which `S_max` — `2` (local hidden variables), `2√2`
(quantum, Tsirelson [Tsirelson1980](@cite)), `4` (no-signalling,
[PopescuRohrlich1994](@cite)) — is a property of the fetched
[`CHSHBound`](@ref), not of this statement: the inequality is the same one in
every regime, which is exactly why the regime belongs on the value.

`S` is untyped: it is a correlator someone measured, and no quantity names that.

Variables: `S` (the bounded one), `S_max`.
"""
@bound :quantum CHSHInequality(S <= S_max::CHSHBound)

"""
    MerminInequality <: AbstractInequality

The Mermin three-party inequality (Mermin, [Mermin1990](@cite)): a measured
Mermin operator value cannot exceed what the theory admits,

`M ≤ M_max`

(slack `M_max − M`), with `M_max = 2` under local realism and `4` in quantum
mechanics — saturated by the GHZ state, whose violation grows exponentially in
the number of parties.  Regime selection lives on [`MerminGHZBound`](@ref), as
for [`CHSHInequality`](@ref).

Variables: `M` (the bounded one), `M_max`.
"""
@bound :quantum MerminInequality(M <= M_max::MerminGHZBound)

# ─── Dynamical bounds: chaos, speed, scrambling ─────────────────────────

"""
    LyapunovChaosBound <: AbstractInequality

The Maldacena–Shenker–Stanford bound on quantum chaos
([MaldacenaShenkerStanford2016](@cite)): the Lyapunov exponent extracted from
the exponential growth of an out-of-time-order correlator cannot exceed the
thermal ceiling,

`λ_L ≤ λ_max = 2π/β`

(slack `λ_max − λ_L`; `ħ = k_B = 1`).  Saturation is the diagnostic of maximal
chaos — holographic duals and large-`N` SYK sit on the bound — but that a given
model saturates it is a model-specific claim and belongs on a consumer's
registry row, not here.

Variables: `λ_L` (the bounded one), `λ_max`.
"""
@bound :quantum LyapunovChaosBound(λ_L <= λ_max::ChaosBound)

"""
    OrthogonalizationTimeBound <: AbstractInequality

The quantum speed limit as a bound on a *measured* evolution time: the time in
which a state actually reaches an orthogonal one cannot be shorter than the
limit its energy data allows,

`τ ≥ τ_min`

(slack `τ − τ_min`), with `τ_min` the fetched [`QuantumSpeedLimit`](@ref) —
the tighter of Margolus–Levitin ([MargolusLevitin1998](@cite)) and
Mandelstam–Tamm.

This is the *value-fetching* form.  [`MargolusLevitinBound`](@ref) and
[`MandelstamTammBound`](@ref) state the same physics directly from `E_above`
and `ΔE`, with no fetched bound — use those when you have the energy data and
this one when you have the limit.

Variables: `τ` (the bounded one), `τ_min`.
"""
@bound :quantum OrthogonalizationTimeBound(τ >= τ_min::QuantumSpeedLimit)

"""
    FastScramblingBound <: AbstractInequality

The fast-scrambling conjecture (Sekino & Susskind, [SekinoSusskind2008](@cite)):
no thermal system scrambles local information into global entanglement faster
than

`t_scr ≥ t_* = (β/2π) log N`

(slack `t_scr − t_*`), with `t_*` the fetched [`ScramblingTime`](@ref).  Black
holes are conjectured to saturate it, which is what makes them the fastest
scramblers in nature.

A conjecture, not a theorem — stated here because it is universal in form and
falsifiable by a measured `t_scr`, which is precisely what a bound in this
package is for.

Variables: `t_scr` (the bounded one), `t_min`.
"""
@bound :quantum FastScramblingBound(t_scr >= t_min::ScramblingTime)

# ─── Quantum-information operational bounds ─────────────────────────────

"""
    SecretKeyRateBound <: AbstractInequality

The BB84 secret-key rate is ACHIEVABLE, so it bounds the extractable key
fraction from BELOW (Shor & Preskill, [ShorPreskill2000](@cite)):

`r ≥ r_min = 1 − 2 H₂(e)`

(slack `r − r_min`), with `r_min` the fetched [`BB84KeyRate`](@ref) at qubit
error rate `e`, positive for `e < 11%`.

The direction is the one that is easy to state backwards: a protocol achieving
*less* than the proven rate is the failure, not one achieving more.

Variables: `r` (the bounded one), `r_min`.
"""
@bound :quantum SecretKeyRateBound(r >= r_min::BB84KeyRate)

"""
    CloningFidelityBound <: AbstractInequality

The no-cloning theorem, quantitatively (Bužek & Hillery,
[BuzekHillery1996](@cite)): no universal `1 → 2` qubit cloner achieves a
single-copy fidelity above the optimal one,

`F ≤ F_max = 5/6`

(slack `F_max − F`), with `F_max` the fetched [`OptimalCloningFidelity`](@ref).
A reported cloner fidelity above it is an error in the calculation, not a
discovery.

Variables: `F` (the bounded one), `F_max`.
"""
@bound :quantum CloningFidelityBound(F <= F_max::OptimalCloningFidelity)

# ─── Holographic ────────────────────────────────────────────────────────

"""
    BekensteinEntropyBound <: AbstractInequality

The Bekenstein universal entropy bound ([Bekenstein1981](@cite)): the entropy of
a system confined to a region of radius `R` with total energy `E` cannot exceed

`S ≤ S_max = 2π R E`

(slack `S_max − S`; `ħ = c = k_B = 1`), with `S_max` the fetched
[`BekensteinBound`](@ref).  Saturated by a black hole, where it reduces to the
area law — the statement that made entropy a geometric quantity.

Variables: `S` (the bounded one), `S_max`.
"""
@bound :holographic BekensteinEntropyBound(S <= S_max::BekensteinBound)
