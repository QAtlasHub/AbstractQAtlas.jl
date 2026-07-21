# relations/entanglement.jl — entanglement-entropy relations.
#
# The identities a many-body calculation (MPS/ED) checks its measured
# entanglement against: the Rényi–purity link, the 1D-CFT logarithmic
# growth that reads off the central charge, and Page's average-entropy
# formula for a random pure state.

"""
    RenyiTwoPurity <: AbstractRelation

The Rényi-2 entanglement entropy as (minus) the log purity,

`S_2 = −ln Tr(ρ_A²) = −ln(purity)`

— the `n = 2` member of `S_n = (1−n)⁻¹ ln Tr ρ_A^n`, the one directly
accessible from the [`Purity`](@ref).

Variables: `S2`, `purity`.
"""
@relation :entanglement RenyiTwoPurity(S2, purity) = S2 + log(purity)

"""
    CFTEntanglementSlope <: AbstractRelation

The logarithmic growth of the entanglement entropy of an interval in a
1D conformal field theory reads off the central charge (Calabrese &
Cardy, [CalabreseCardy2004](@cite)): for a subsystem of size `ℓ`,

`S(ℓ) = (c/3) ln ℓ + const`   ⟹   `dS/d(ln ℓ) = c/3`   (periodic BC;
the open-boundary coefficient is `c/6`).

Supplied-derivative convention: `dS_dlogℓ` is the caller-computed slope
of `S` against `ln ℓ`.  Variables: `dS_dlogℓ`, `c`.
"""
@relation :entanglement CFTEntanglementSlope(dS_dlogℓ, c) = dS_dlogℓ - c / 3

"""
    page_average_entropy(dA, dB) -> Float64

Page's average entanglement entropy of the smaller subsystem `A` for a
Haar-random pure state of a bipartite system `A ⊗ B` with Hilbert-space
dimensions `dA ≤ dB` (Page, [Page1993](@cite)):

`⟨S_A⟩ = ( Σ_{k=dB+1}^{dA·dB} 1/k ) − (dA − 1)/(2 dB)`.

Nearly maximal, `⟨S_A⟩ ≈ ln dA − dA/(2 dB)`: a random state is almost
maximally entangled, deficit `dA/(2dB)`.  `dA > dB` is symmetric — call
with the smaller dimension first.
"""
function page_average_entropy(dA::Integer, dB::Integer)
    (dA >= 1 && dB >= 1) || error("dimensions must be ≥ 1")
    dA <= dB || return page_average_entropy(dB, dA)   # symmetric; use the smaller
    harmonic = sum(1 / k for k in (dB + 1):(dA * dB))
    return harmonic - (dA - 1) / (2 * dB)
end
export page_average_entropy

# ─── Entropy / quantum-information inequalities (≥ 0 slack; @inequality) ──
#
# The bound-type constraints a many-body entanglement calculation must
# satisfy — the first users of the AbstractInequality kind.  Each holds
# iff its slack `≥ 0`; `check` tests that direction, `solve` gives the
# saturation (tight-bound) value.

"""
    EntropyNonNegativity <: AbstractInequality

The von Neumann / Rényi entanglement entropy is non-negative, `S ≥ 0`
(slack `S`).  Saturated by a pure (unentangled) subsystem.

Variables: `S`.
"""
@inequality :entanglement EntropyNonNegativity(S) = S

"""
    MaxEntropyBound <: AbstractInequality

The entropy of a subsystem cannot exceed the log of its Hilbert-space
dimension, `S ≤ ln d` (slack `ln d − S`).  Saturated by the maximally
mixed state; the gap `ln d − S` is the maximal-entanglement deficit.

Variables: `S`, `log_d` = `ln d`.
"""
@inequality :entanglement MaxEntropyBound(S, log_d) = log_d - S

"""
    Subadditivity <: AbstractInequality

Subadditivity of the von Neumann entropy, `S(AB) ≤ S(A) + S(B)` (slack
`S_A + S_B − S_AB` — the mutual information `I(A:B) ≥ 0`; Araki & Lieb,
[ArakiLieb1970](@cite)).  Saturated by a product state
`ρ_AB = ρ_A ⊗ ρ_B`.

Variables: `S_A`, `S_B`, `S_AB`.
"""
@inequality :entanglement Subadditivity(S_A, S_B, S_AB) = S_A + S_B - S_AB

"""
    ArakiLieb <: AbstractInequality

The Araki–Lieb triangle inequality, `S(AB) ≥ |S(A) − S(B)|` (slack
`S_AB − |S_A − S_B|`; Araki & Lieb, [ArakiLieb1970](@cite)) —
the lower companion of [`Subadditivity`](@ref).  Saturated when one
subsystem purifies the other.

Variables: `S_AB`, `S_A`, `S_B`.
"""
@inequality :entanglement ArakiLieb(S_AB, S_A, S_B) = S_AB - abs(S_A - S_B)

"""
    StrongSubadditivity <: AbstractInequality

Strong subadditivity of the quantum entropy,
`S(ABC) + S(B) ≤ S(AB) + S(BC)` (slack `S_AB + S_BC − S_ABC − S_B`; Lieb
& Ruskai, [LiebRuskai1973](@cite)) — equivalently the conditional
mutual information `I(A:C|B) ≥ 0`.  The deepest entropy inequality; the
monogamy backbone of quantum information.

Variables: `S_AB`, `S_BC`, `S_ABC`, `S_B`.
"""
@inequality :entanglement StrongSubadditivity(S_AB, S_BC, S_ABC, S_B) =
    S_AB + S_BC - S_ABC - S_B

"""
    WeakMonotonicity <: AbstractInequality

Weak monotonicity of the quantum entropy,
`S(A) + S(C) ≤ S(AB) + S(BC)` (slack `S_AB + S_BC − S_A − S_C`) — the
purification dual of [`StrongSubadditivity`](@ref) (purify `C`; SSA on the
purified state *is* weak monotonicity here), equivalent and equally universal,
but stated in the *outer* regions `A, C` rather than `ABC, B`.  Requires strictly
less than SSA — no full-system `S(ABC)` — so it is checkable from partial data.

Variables: `S_AB`, `S_BC`, `S_A`, `S_C`.
"""
@inequality :entanglement WeakMonotonicity(S_AB, S_BC, S_A, S_C) = S_AB + S_BC - S_A - S_C

"""
    RenyiMonotonicity <: AbstractInequality

The Rényi entropy `S_α` is non-increasing in the order `α`: for
`α_low < α_high`, `S_{α_low} ≥ S_{α_high}` (slack `S_low − S_high`).  In
particular `S_0 ≥ S_1 (von Neumann) ≥ S_2 ≥ … ≥ S_∞`.

Variables: `S_low` = `S_{α_low}`, `S_high` = `S_{α_high}` (with `α_low < α_high`).
"""
@inequality :entanglement RenyiMonotonicity(S_low, S_high) = S_low - S_high

# ─── The entropy zoo: Rényi / Tsallis / mutual / conditional / relative ──
#
# The defining relations of the one-parameter entropy families and the
# composite (multi-party) entropies — unifying the entanglement measures a
# many-body calculation reports (issue #27).

"""
    RenyiEntropyMoment <: AbstractRelation

The Rényi entropy from the density-matrix moment `Tr ρ^α` (α ≠ 1),

`S_α = ln(Tr ρ^α) / (1 − α)`,

the general form behind [`RenyiTwoPurity`](@ref) (α = 2: `S_2 = −ln Tr ρ²`)
and, as `α → 1`, the [`VonNeumannEntropy`](@ref).

Variables: `Sα`, `moment` = `Tr ρ^α`, `α`.
"""
@relation :entanglement RenyiEntropyMoment(Sα, moment, α) = Sα - log(moment) / (1 - α)

"""
    TsallisEntropyMoment <: AbstractRelation

The Tsallis entropy from the moment `Tr ρ^q` (q ≠ 1; Tsallis, [Tsallis1988](@cite)),

`S_q = (1 − Tr ρ^q) / (q − 1)`.

Variables: `Sq`, `moment` = `Tr ρ^q`, `q`.
"""
@relation :entanglement TsallisEntropyMoment(Sq, moment, q) = Sq - (1 - moment) / (q - 1)

"""
    MutualInformationDefinition <: AbstractRelation

The quantum mutual information,

`I(A:B) = S(A) + S(B) − S(AB)`,

([`MutualInformation`](@ref); non-negative by [`Subadditivity`](@ref)).

Variables: `I`, `S_A`, `S_B`, `S_AB`.
"""
@relation :entanglement MutualInformationDefinition(I, S_A, S_B, S_AB) =
    I - (S_A + S_B - S_AB)

"""
    ConditionalEntropyDefinition <: AbstractRelation

The quantum conditional entropy,

`S(A|B) = S(AB) − S(B)`,

([`ConditionalEntropy`](@ref); can be negative — an entanglement witness).

Variables: `S_cond`, `S_AB`, `S_B`.
"""
@relation :entanglement ConditionalEntropyDefinition(S_cond, S_AB, S_B) =
    S_cond - (S_AB - S_B)

"""
    RelativeEntropyNonNegativity <: AbstractInequality

Klein's inequality: the quantum relative entropy is non-negative,

`S(ρ‖σ) ≥ 0`,

(slack `S_rel`; zero iff `ρ = σ`).  The bedrock positivity behind
subadditivity and the second law (Lindblad, [Lindblad1975](@cite); Vedral, [Vedral2002](@cite)).

Variables: `S_rel` = `S(ρ‖σ)`.
"""
@inequality :entanglement RelativeEntropyNonNegativity(S_rel) = S_rel

# ─── Entropy of mixing: concavity + the Holevo upper bound ───────────────
#
# For a mixture ρ = Σᵢ pᵢ ρᵢ the entropy is sandwiched by the weighted-average
# component entropy Σᵢ pᵢ S(ρᵢ): concavity from below, and from above by that
# average plus the classical mixing entropy H(p).  The gap S(ρ) − Σᵢ pᵢ S(ρᵢ)
# is the Holevo χ, `0 ≤ χ ≤ H(p)`.  `S_avg` and `H_weights` are caller-supplied
# aggregates over the ensemble (the sum over members is a functional step).

"""
    EntropyMixingConcavity <: AbstractInequality

Concavity of the von Neumann entropy — mixing states cannot decrease the entropy,

`S(Σᵢ pᵢ ρᵢ) ≥ Σᵢ pᵢ S(ρᵢ)`

(slack `S_mix − S_avg`; Wehrl, [Wehrl1978](@cite)).  Saturated when every
`ρᵢ` with `pᵢ > 0` is the same state.

Variables: `S_mix` = `S(Σᵢ pᵢ ρᵢ)`, `S_avg` = the caller-supplied `Σᵢ pᵢ S(ρᵢ)`.
"""
@inequality :entanglement EntropyMixingConcavity(S_mix, S_avg) = S_mix - S_avg

"""
    HolevoMixingBound <: AbstractInequality

The upper companion of [`EntropyMixingConcavity`](@ref): the entropy of a mixture
exceeds the average component entropy by at most the classical mixing entropy,

`S(Σᵢ pᵢ ρᵢ) ≤ Σᵢ pᵢ S(ρᵢ) + H(p)`,   `H(p) = −Σᵢ pᵢ ln pᵢ`

(slack `S_avg + H_weights − S_mix`; Wehrl, [Wehrl1978](@cite)).  Saturated
when the `ρᵢ` have mutually orthogonal support; the gap `S_mix − S_avg` is the Holevo
`χ`, bounded in `[0, H(p)]`.

Variables: `S_avg` = `Σᵢ pᵢ S(ρᵢ)`, `H_weights` = `H(p)`, `S_mix` = `S(Σᵢ pᵢ ρᵢ)`.
"""
@inequality :entanglement HolevoMixingBound(S_avg, H_weights, S_mix) =
    S_avg + H_weights - S_mix

# ─── Measurement and quantum-Markov entropies ───────────────────────────

"""
    MeasurementEntropyIncrease <: AbstractInequality

A projective measurement (dephasing) does not decrease the entropy,

`S(Δρ) ≥ S(ρ)`

(slack `S_meas − S`; [`MeasurementEntropy`](@ref)).  Saturated iff `ρ` is
already diagonal in the measurement basis (`Δρ = ρ`).

Variables: `S_meas` = `S(Δρ)`, `S` = `S(ρ)`.
"""
@inequality :entanglement MeasurementEntropyIncrease(S_meas, S) = S_meas - S

"""
    MeasurementEntropyRelative <: AbstractRelation

The entropy gain from a projective measurement equals the relative entropy
to the dephased state,

`S(Δρ) − S(ρ) = S(ρ‖Δρ)`,

tying the [`MeasurementEntropy`](@ref) to the [`RelativeEntropy`](@ref)
(Vedral, [Vedral2002](@cite)).

Variables: `S_meas` = `S(Δρ)`, `S` = `S(ρ)`, `S_rel` = `S(ρ‖Δρ)`.
"""
@relation :entanglement MeasurementEntropyRelative(S_meas, S, S_rel) = (S_meas - S) - S_rel

"""
    MarkovEntropyDefinition <: AbstractRelation

The conditional mutual information (the [`MarkovEntropy`](@ref)),

`I(A:C|B) = S(AB) + S(BC) − S(ABC) − S(B)`,

equal to the strong-subadditivity slack ([`StrongSubadditivity`](@ref));
its vanishing marks a quantum Markov chain `A–B–C` (Hayden, Jozsa, Petz &
Winter, [HaydenJozsaPetzWinter2004](@cite)).

Variables: `I_cmi`, `S_AB`, `S_BC`, `S_ABC`, `S_B`.
"""
@relation :entanglement MarkovEntropyDefinition(I_cmi, S_AB, S_BC, S_ABC, S_B) =
    I_cmi - (S_AB + S_BC - S_ABC - S_B)

# ─── Free-fermion (Gaussian) entanglement from the correlation matrix ────
#
# For a Gaussian (free-fermion) state the reduced density matrix ρ_A is
# fixed ENTIRELY by the restricted two-point correlation matrix
# C_ij = ⟨c†_i c_j⟩|_A — because Wick's theorem (relations/wick.jl) makes
# every higher moment a determinant of C.  ρ_A = e^{−H_ent}/Z with a
# QUADRATIC entanglement Hamiltonian, so the eigenvalues ζ_k ∈ [0,1] of
# C_A give the whole entanglement spectrum (Peschel, [Peschel2003](@cite);
# Chung & Peschel, [ChungPeschel2001](@cite)).  This mapping
# holds for Gaussian states ONLY — an interacting ρ_A is not fixed by its
# two-point function.

"""
    EntanglementSpectrumCorrelation <: AbstractRelation

The free-fermion entanglement (single-particle) spectrum from the
correlation-matrix eigenvalue `ζ ∈ (0, 1)` (Peschel, [Peschel2003](@cite)),

`ε = ln((1 − ζ)/ζ)`,

the eigenvalue of the quadratic entanglement Hamiltonian `H_ent`; inverting
gives the Fermi-Dirac occupation `ζ = 1/(e^ε + 1)`.  A maximally-entangled
mode `ζ = ½` sits at `ε = 0`.  **Gaussian states only** — the correlation
matrix fixes `ρ_A` via Wick's theorem ([`wick_contraction`](@ref)).

Variables: `ε`, `ζ`.
"""
@relation :entanglement EntanglementSpectrumCorrelation(ε, ζ) = ε - log((1 - ζ) / ζ)

"""
    free_fermion_entanglement_entropy(ζ) -> Float64

The von Neumann entanglement entropy of a free-fermion (Gaussian) region
from the eigenvalues `ζ_k ∈ [0, 1]` of its restricted correlation matrix
`C_ij = ⟨c†_i c_j⟩` (Peschel, [Peschel2003](@cite)),

`S_A = −Σ_k [ζ_k ln ζ_k + (1 − ζ_k) ln(1 − ζ_k)]`

(the sum of per-mode binary entropies).  A fully occupied/empty mode
(`ζ = 0, 1`) contributes nothing; a maximally-entangled mode (`ζ = ½`)
contributes `ln 2`.  Valid for **Gaussian states only**.
"""
function free_fermion_entanglement_entropy(ζ)
    h(x) = (x <= 0 || x >= 1) ? 0.0 : -x * log(x) - (1 - x) * log(1 - x)
    return sum(h, ζ)
end
export free_fermion_entanglement_entropy

"""
    free_fermion_renyi_entropy(ζ, n) -> Float64

The order-`n` Rényi entanglement entropy of a free-fermion region from the
correlation-matrix eigenvalues `ζ_k` (`n ≠ 1`),

`S_A^{(n)} = (1 − n)⁻¹ Σ_k ln[ζ_k^n + (1 − ζ_k)^n]`,

recovering [`free_fermion_entanglement_entropy`](@ref) as `n → 1`.
"""
function free_fermion_renyi_entropy(ζ, n)
    n == 1 && error(
        "Rényi order n = 1 is the von Neumann limit; use free_fermion_entanglement_entropy",
    )
    return sum(x -> log(x^n + (1 - x)^n), ζ) / (1 - n)
end
export free_fermion_renyi_entropy

# ─── Multipartite entanglement: monogamy, tangle, tripartite (#27) ───────

"""
    ConcurrenceTangle <: AbstractRelation

The tangle is the squared concurrence (Wootters, [Wootters1998](@cite)),

`τ = C²`,

([`Tangle`](@ref), [`Concurrence`](@ref)).

Variables: `τ`, `C`.
"""
@relation :entanglement ConcurrenceTangle(τ, C) = τ - C^2

"""
    Monogamy <: AbstractInequality

The Coffman–Kundu–Wootters monogamy of entanglement (Coffman, Kundu &
Wootters, [CoffmanKunduWootters2000](@cite)): the tangle of `A` with the rest
bounds the sum of its pairwise tangles,

`τ(A:BC) ≥ τ(A:B) + τ(A:C)`

(slack `τ_ABC − τ_AB − τ_AC` = the [`ThreeTangle`](@ref) `τ₃ ≥ 0`).
Entanglement cannot be freely shared.

Variables: `τ_ABC`, `τ_AB`, `τ_AC`.
"""
@inequality :entanglement Monogamy(τ_ABC, τ_AB, τ_AC) = τ_ABC - τ_AB - τ_AC

"""
    ThreeTangleDefinition <: AbstractRelation

The residual three-tangle — the genuinely tripartite entanglement beyond
the pairwise budget (Coffman, Kundu & Wootters, [CoffmanKunduWootters2000](@cite)),

`τ₃ = τ(A:BC) − τ(A:B) − τ(A:C)`,

the saturation gap of [`Monogamy`](@ref).

Variables: `τ3`, `τ_ABC`, `τ_AB`, `τ_AC`.
"""
@relation :entanglement ThreeTangleDefinition(τ3, τ_ABC, τ_AB, τ_AC) =
    τ3 - (τ_ABC - τ_AB - τ_AC)

"""
    TripartiteInformationDefinition <: AbstractRelation

The tripartite information,

`I₃(A:B:C) = I(A:B) + I(A:C) − I(A:BC)`,

([`TripartiteInformation`](@ref)); a negative `I₃` signals genuinely
multipartite (scrambled) correlation.

Variables: `I3`, `I_AB`, `I_AC`, `I_ABC`.
"""
@relation :entanglement TripartiteInformationDefinition(I3, I_AB, I_AC, I_ABC) =
    I3 - (I_AB + I_AC - I_ABC)

"""
    KitaevPreskillTEE <: AbstractRelation

The topological entanglement entropy from a tripartition (Kitaev &
Preskill, [KitaevPreskill2006](@cite)),

`S_A + S_B + S_C − S_AB − S_BC − S_CA + S_ABC = −γ`,

the universal constant `γ = ln D` isolated from the area law by the
alternating tripartite sum (`γ > 0` ⇒ topological order).

Variables: `γ`, `S_A`, `S_B`, `S_C`, `S_AB`, `S_BC`, `S_CA`, `S_ABC`.
"""
@relation :entanglement KitaevPreskillTEE(γ, S_A, S_B, S_C, S_AB, S_BC, S_CA, S_ABC) =
    (S_A + S_B + S_C - S_AB - S_BC - S_CA + S_ABC) + γ
