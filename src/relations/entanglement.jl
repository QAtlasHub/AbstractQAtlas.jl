# relations/entanglement.jl — entanglement-entropy relations.
#
# The identities a many-body calculation (MPS/ED) checks its measured
# entanglement against: the Rényi–purity link, the 1D-CFT logarithmic
# growth that reads off the central charge, and Page's average-entropy
# formula for a random pure state.
#
# Entropies here are in NATS. The entanglement literature usually counts bits,
# `S = -Tr ρ log₂ ρ`, and the two differ by `ln 2`. A coefficient multiplying a
# logarithm (`c`, `c̃`) is unaffected, since the base rescales entropy and
# logarithm alike; an additive constant or a bare difference of entropies is not,
# and that is where a transcribed formula silently gains or loses a factor.

# The conformal chord the finite-size forms are stated in terms of. Shared so the
# relations and the Region-typed layer above them cannot drift apart.
_chord(L, ℓ) = (L / π) * sin(π * ℓ / L)

# `sin` is periodic, so an ℓ outside the chain does not merely give a wrong number:
# ℓ = 250 on L = 100 returns exactly the ℓ = 50 answer. And ℓ = L is not caught by a
# blow-up either, since `sin(float(π))` is 1.2e-16 rather than 0, which turns into a
# large finite value dominated by rounding. Both are refused.
function _require_block(what::Symbol, L, ℓ)
    return 0 < ℓ < L || error(
        "$what: need 0 < ℓ < L, got ℓ = $ℓ on L = $L. Outside that range `sin(πℓ/L)` " *
        "aliases onto a legitimate answer; at ℓ = L it is rounding noise, where the " *
        "true entropy of the whole system is 0.",
    )
end

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

Logarithmic growth of a region's entanglement entropy in a 1D CFT, reading off
the central charge (Calabrese & Cardy, [CalabreseCardy2004](@cite)):

`dS/d(ln ℓ) = ncuts · c/6`.

`ncuts` counts the cuts bounding the region — set by where it sits, not by the
chain's boundary condition:

| region | `ncuts` | |
|---|---|---|
| one interval on a ring | 2 | `c/3` |
| block at an open end | 1 | `c/6` |
| block in the bulk of an open chain | 2 | `c/3` |

The last row is the one "OBC ⇒ `c/6`" gets wrong.

Variables: `dS_dlogℓ` (caller-computed slope against `ln ℓ`), `c`, `ncuts`.
`c` is the typed subject; [`VonNeumannEntropy`](@ref) arrives via the supplied
derivative, hence [`also_constrains`](@ref).  `ncuts` has no default:
[`Region`](@ref) carries no adjacency or boundary, so nothing can compute it.
"""
@relation :entanglement CFTEntanglementSlope(dS_dlogℓ, c::CentralCharge, ncuts) =
    dS_dlogℓ - ncuts * c / 6

"""
    InfiniteRandomnessEntanglementSlope <: AbstractRelation

The same logarithmic growth at a **one-dimensional** infinite-randomness fixed
point, with the effective central charge in place of the CFT one (Refael &
Moore, [RefaelMoore2004](@cite)):

`dS/d(ln ℓ) = ncuts · c̃/6`.

No `d` slot, because there is no family to index: above one dimension the entropy
obeys an area law, and whether the fixed point is reached at all is
model-dependent and settled only numerically, the random transverse-field Ising
model reaching one in `d ≥ 2` where the random Heisenberg antiferromagnet does not
([IgloiMonthus2005](@cite), Sec. 9).

The claim is the leading slope, not the finite-size form: see
[`CFTEntanglementPBC`](@ref) and [`CFTEntanglementOBC`](@ref) for what else a
boundary changes, none of which follows from a slope at a fixed point with no
conformal map.

The fixed point is not conformally invariant, yet both the form and the geometry
factor survive: `ncuts` means what it means in [`CFTEntanglementSlope`](@ref),
and only `c` becomes [`EffectiveCentralCharge`](@ref).  Refael and Moore's `2 ln √L` (Eq. 13) is that same
factor arriving from the RG, two cuts each contributing at `Γ = √L`; the source
corrects that equation's rate three sentences later, and it is the two, which
survives into Eq. (19), that is being leaned on here.

`c̃` is measured per class, not derived: `(ln 2)/2` for the random transverse
field Ising chain, which grows as `(ln 2/6) ln ℓ ≈ 0.1155 ln ℓ` across two cuts,
and `ln 2` for the random singlet phase of the Heisenberg and XX chains.  Both
are `ln 2` times the pure value, which the source reports for every chain it
treats while calling a general law only possible, so that is not a relation
here.  `c̃` is the same number in either base (see the file header).

Variables: `dS_dlogℓ`, `c̃`, `ncuts`; the entropy arrives through the supplied
derivative, hence [`also_constrains`](@ref).  With a region and a boundary
condition in hand, derive `ncuts` from [`entanglement_cuts`](@ref) rather than
passing a literal.
"""
@relation :entanglement InfiniteRandomnessEntanglementSlope(
    dS_dlogℓ, c̃::EffectiveCentralCharge, ncuts
) = dS_dlogℓ - ncuts * c̃ / 6

"""
    entanglement_cuts(bc::BoundaryCondition, A::Region) -> Int

The number of cuts bounding `A`, which Iglói & Lin ([IgloiLin2008](@cite),
Eq. 5) call `b`, "the number of boundary points between the subsystem and the
rest of the chain": the count of adjacent site pairs with exactly one member in
`A`, with `(N, 1)` adjacent under [`PBC`](@ref).

This is the `ncuts` the entanglement relations take, derived instead of
asserted; [`CFTEntanglementSlope`](@ref)'s table says what hand-supplying it
gets wrong.

[`Region`](@ref) is a set with no adjacency, so the sites must be integers and
the chain length must come from `bc`; anything else is refused rather than
guessed.  A block filling the whole ring returns 0, as it must.

```julia
entanglement_cuts(PBC(8), Region(2, 3, 4))   # 2
entanglement_cuts(OBC(8), Region(1, 2, 3))   # 1, it touches the end
entanglement_cuts(OBC(8), Region(2, 3, 4))   # 2, the same block in the bulk
```
"""
function entanglement_cuts(bc::BoundaryCondition, A::Region{<:Integer})
    isempty(A) && return 0
    eltype(A.sites) === Bool && error(
        "entanglement_cuts: Bool is an Integer but not a lattice index; got $(A.sites)."
    )
    # Signed arithmetic: `lo - 1` on unsigned labels wraps to typemax and makes the
    # range below empty, which would return 0 cuts for a region that has two.
    sites = Set{Int}(Int(i) for i in A.sites)
    if bc isa Infinite
        lo, hi = minimum(sites), maximum(sites)
        return count(i -> (i in sites) != (i + 1 in sites), (lo - 1):hi)
    end
    bc isa Union{OBC,PBC} || error(
        "entanglement_cuts: no adjacency defined for $(typeof(bc)). Add a method here " *
        "when a boundary condition is added, rather than letting it read as open.",
    )
    N = bc.N
    N > 0 || error(
        "entanglement_cuts: $bc declares no chain length. Pass it as OBC(N) / PBC(N); " *
        "the `N = 0` sentinel means the size lives in a caller's kwargs, which this " *
        "function cannot see.",
    )
    maximum(sites) <= N && minimum(sites) >= 1 || error(
        "entanglement_cuts: sites $(minimum(sites))..$(maximum(sites)) fall outside " *
        "the chain 1..$N declared by $bc.",
    )
    n = count(i -> (i in sites) != (i + 1 in sites), 1:(N - 1))
    bc isa PBC && ((N in sites) != (1 in sites)) && (n += 1)
    return n
end

# `Region` is a set layer with no adjacency, so a non-integer label has no
# neighbour to be separated from. Refused here rather than by a MethodError, which
# would not say what to do instead.
function entanglement_cuts(::BoundaryCondition, A::Region)
    isempty(A) && return 0
    return error(
        "entanglement_cuts: adjacency needs integer sites; got $(eltype(A.sites)). " *
        "Region is a set layer with no geometry, so pass `ncuts` directly instead.",
    )
end
export entanglement_cuts

"""
    CFTEntanglementPBC <: AbstractRelation

Entanglement entropy of a block of `ℓ` sites in a critical ring of `L`
(Iglói & Lin, [IgloiLin2008](@cite), Eq. 2):

`S = (c/3) ln[(L/π) sin(πℓ/L)] + c₁`.

Two cuts, hence `c/3`.  `c₁` is not universal and moves with the base (see the
file header).  As `ℓ ≪ L` the chord tends to `ℓ` and this becomes the
infinite-chain `S = (c/3) ln ℓ + c₁` of Eq. (4).
"""
@relation :entanglement CFTEntanglementPBC(S, c::CentralCharge, L, ℓ, c₁) = begin
    _require_block(:CFTEntanglementPBC, L, ℓ)
    S - (c / 3) * log(_chord(L, ℓ)) - c₁
end

"""
    CFTEntanglementInfinite <: AbstractRelation

The thermodynamic limit of [`CFTEntanglementPBC`](@ref) (Iglói & Lin,
[IgloiLin2008](@cite), Eq. 4):

`S = (c/3) ln ℓ + c₁`.

Still two cuts, and the same `c₁`: the chord tends to `ℓ` as `ℓ ≪ L`, so this is
where the ring form goes rather than a separate law.  It exists as its own
relation because a bag on an infinite chain has no `L` to supply.
"""
@relation :entanglement CFTEntanglementInfinite(S, c::CentralCharge, ℓ, c₁) = begin
    ℓ > 0 || error("CFTEntanglementInfinite: need ℓ > 0, got $ℓ.")
    S - (c / 3) * log(ℓ) - c₁
end

"""
    CFTEntanglementOBC <: AbstractRelation

The same for the leftmost `ℓ` sites of a critical open chain of `L`
(Iglói & Lin, [IgloiLin2008](@cite), Eq. 3):

`S = (c/6) ln[(2L/π) sin(πℓ/L)] + ln g + c₁/2`.

Three things separate this from [`CFTEntanglementPBC`](@ref), and only the
first is the cut count: one cut gives `c/6`, the chord carries `2L/π` rather
than `L/π`, and an open chain has a boundary entropy `ln g` (Affleck & Ludwig)
that a ring does not.  `c₁` is the same constant as in the ring, entering
halved, so reading one geometry's data with the other's formula misses all
three.
"""
@relation :entanglement CFTEntanglementOBC(S, c::CentralCharge, L, ℓ, c₁, ln_g) = begin
    _require_block(:CFTEntanglementOBC, L, ℓ)
    S - (c / 6) * log(2 * _chord(L, ℓ)) - ln_g - c₁ / 2
end

"""
    OffCriticalEntanglementSaturation <: AbstractRelation

Away from criticality the entropy stops growing with `ℓ` and saturates on the
correlation length (Iglói & Lin, [IgloiLin2008](@cite), Eq. 5, valid for
`ξ ≪ ℓ`):

`S ≃ ncuts · (c/6) ln ξ`.

The source's `b` is `ncuts` (see [`entanglement_cuts`](@ref)): the same axis
counts here, with `ξ` in the place `ℓ` held.
"""
@relation :entanglement OffCriticalEntanglementSaturation(
    S, c::CentralCharge, ξ::CorrelationLength, ncuts
) = S - ncuts * (c / 6) * log(ξ)

"""
    HalvedChainEntropyDifference <: AbstractRelation

Central charge from two chain lengths rather than from a fit (Iglói & Lin,
[IgloiLin2008](@cite), Sec. 3.1), with `ΔS(L) = S_L(L/2) - S_{L/2}(L/4)`:

`ΔS = ncuts · (c/6) ln 2`.

Exact on the conformal forms, since halving the chain shifts the chord by a
factor of two and the non-universal `c₁` cancels: the estimator needs no
constant, which is why the source uses it.  The source reads `ΔS = c/3` for a
ring and `c/6` for an open chain because it counts bits, where `log₂ 2 = 1`
absorbs the factor; in nats it does not, and dropping it returns `c ln 2` for
`c`.  For the Ising chain that is `(ln 2)/2`, which is
exactly the effective central charge of the *random* Ising chain, so the slip is
numerically indistinguishable from having measured a different fixed point.

The finite-size approach differs by boundary condition in its exponent, not only
its amplitude: the source measures `c(L) = 1/2 - 0.623/L² + O(L⁻³)` on a ring
against `c(L) = 1/2 + 1.339/L + O(L⁻²)` on an open chain, so an open chain's
leading correction is one power of `L` slower.
"""
@relation :entanglement HalvedChainEntropyDifference(ΔS, c::CentralCharge, ncuts) =
    ΔS - ncuts * (c / 6) * log(2)

"""
    InfiniteRandomnessEntanglementPBC <: AbstractRelation

Finite-size entropy of a block in a random critical ring, at the
infinite-randomness fixed point (Iglói & Lin, [IgloiLin2008](@cite), Eq. 24):

`S̄ = (c̃/3) ln[L f(ℓ/L)] + c₁′`.

`f` is caller-supplied and is not the conformal chord.  It is reflection
symmetric, tends to `v` as `v → 0`, and expands as `f(v) = Σₖ Aₖ sin((2k-1)πv)`
under `Σₖ Aₖ(2k-1)π = 1`, where the source notes that a conformally invariant
model has **only the first term**.  Keeping `k = 1` forces `A₁ = 1/π` and returns
`L f = (L/π) sin(πℓ/L)`, which is [`CFTEntanglementPBC`](@ref) exactly, so at
finite size the two differ by the higher harmonics and not by the coefficient
alone.

Only the value of `f` reaches the relation, so `L > 0` and `f > 0` are all it can
check and a small positive `f` still returns a large negative entropy; reached
through [`finite_size_entropy_report`](@ref), `f` is a callable sampled at
`v = ℓ/L` from the region, which is where that is visible.

`c̃ = (ln 2)/2` is reported universal here in a stronger sense than the slope
alone requires, being independent of the form of the disorder, while `c₁′`
depends on it.
"""
@relation :entanglement InfiniteRandomnessEntanglementPBC(
    S̄, c̃::EffectiveCentralCharge, L, f, c₁′
) = begin
    (L > 0 && f > 0) || error(
        "InfiniteRandomnessEntanglementPBC: L = $L and f = $f must both be positive. " *
        "`log(L*f)` is finite whenever their product is, so a negative pair returns " *
        "an ordinary-looking number, and an f near zero returns a negative entropy.",
    )
    S̄ - (c̃ / 3) * log(L * f) - c₁′
end

@experimental """
Eq. (24)'s scaling function is not settled here: `f` is caller-supplied, no
independent oracle checks the form, and the only case pinned is its conformal
one-harmonic reduction. No bound on `f` follows from the geometry either, so a
positive but small value returns a negative entropy rather than a refusal
""" InfiniteRandomnessEntanglementPBC

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

# ─── Entropy / quantum-information bounds (≥ 0 slack; @bound) ───────────
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
@bound :entanglement EntropyNonNegativity(S >= 0)

"""
    MaxEntropyBound <: AbstractInequality

The entropy of a subsystem cannot exceed the log of its Hilbert-space
dimension, `S ≤ ln d` (slack `ln d − S`).  Saturated by the maximally
mixed state; the gap `ln d − S` is the maximal-entanglement deficit.

Variables: `S`, `log_d` = `ln d`.
"""
@bound :entanglement MaxEntropyBound(S <= log_d)

"""
    Subadditivity <: AbstractInequality

Subadditivity of the von Neumann entropy, `S(AB) ≤ S(A) + S(B)` (slack
`S_A + S_B − S_AB` — the mutual information `I(A:B) ≥ 0`; Araki & Lieb,
[ArakiLieb1970](@cite)).  Saturated by a product state
`ρ_AB = ρ_A ⊗ ρ_B`.

Variables: `S_A`, `S_B`, `S_AB`.
"""
@bound :entanglement Subadditivity(S_A, S_B, S_AB) = S_AB <= S_A + S_B

"""
    ArakiLieb <: AbstractInequality

The Araki–Lieb triangle inequality, `S(AB) ≥ |S(A) − S(B)|` (slack
`S_AB − |S_A − S_B|`; Araki & Lieb, [ArakiLieb1970](@cite)) —
the lower companion of [`Subadditivity`](@ref).  Saturated when one
subsystem purifies the other.

Variables: `S_AB`, `S_A`, `S_B`.
"""
@bound :entanglement ArakiLieb(S_AB, S_A, S_B) = S_AB >= abs(S_A - S_B)

"""
    StrongSubadditivity <: AbstractInequality

Strong subadditivity of the quantum entropy,
`S(ABC) + S(B) ≤ S(AB) + S(BC)` (slack `S_AB + S_BC − S_ABC − S_B`; Lieb
& Ruskai, [LiebRuskai1973](@cite)) — equivalently the conditional
mutual information `I(A:C|B) ≥ 0`.  The deepest entropy inequality; the
monogamy backbone of quantum information.

Variables: `S_AB`, `S_BC`, `S_ABC`, `S_B`.
"""
@bound :entanglement StrongSubadditivity(S_AB, S_BC, S_ABC, S_B) =
    S_ABC + S_B <= S_AB + S_BC

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
@bound :entanglement WeakMonotonicity(S_AB, S_BC, S_A, S_C) = S_A + S_C <= S_AB + S_BC

"""
    RenyiMonotonicity <: AbstractInequality

The Rényi entropy `S_α` is non-increasing in the order `α`: for
`α_low < α_high`, `S_{α_low} ≥ S_{α_high}` (slack `S_low − S_high`).  In
particular `S_0 ≥ S_1 (von Neumann) ≥ S_2 ≥ … ≥ S_∞`.

Variables: `S_high` = `S_{α_high}` (the bounded one), `S_low` = `S_{α_low}`
(with `α_low < α_high`).
"""
@bound :entanglement RenyiMonotonicity(S_high <= S_low)

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
@bound :entanglement RelativeEntropyNonNegativity(S_rel >= 0)

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

Variables: `S_avg` = the caller-supplied `Σᵢ pᵢ S(ρᵢ)` (the bounded one),
`S_mix` = `S(Σᵢ pᵢ ρᵢ)`.
"""
@bound :entanglement EntropyMixingConcavity(S_avg <= S_mix)

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
@bound :entanglement HolevoMixingBound(S_avg, H_weights, S_mix) = S_mix <= S_avg + H_weights

# ─── Measurement and quantum-Markov entropies ───────────────────────────

"""
    MeasurementEntropyIncrease <: AbstractInequality

A projective measurement (dephasing) does not decrease the entropy,

`S(Δρ) ≥ S(ρ)`

(slack `S_meas − S`; [`MeasurementEntropy`](@ref)).  Saturated iff `ρ` is
already diagonal in the measurement basis (`Δρ = ρ`).

Variables: `S` = `S(ρ)` (the bounded one), `S_meas` = `S(Δρ)`.
"""
@bound :entanglement MeasurementEntropyIncrease(S <= S_meas)

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
@bound :entanglement Monogamy(τ_ABC, τ_AB, τ_AC) = τ_ABC >= τ_AB + τ_AC

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
