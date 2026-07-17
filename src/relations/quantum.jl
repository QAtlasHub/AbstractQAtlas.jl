# relations/quantum.jl — exact quantum-mechanical identities.
#
# The state-independent identities any quantum (many-body) calculation
# checks its expectation values against: the virial theorem, the
# Hellmann–Feynman and Ehrenfest theorems, the zero-variance eigenstate
# condition (the variational / DMRG convergence metric), the Robertson
# uncertainty inequality, and the dynamical speed limits (Lieb–Robinson
# velocity, Mandelstam–Tamm / Margolus–Levitin orthogonalization time).
# Domain tag :quantum throughout.
#
# References (doiget-verified, docs/references.bib): Feynman, Phys. Rev.
# 56, 340 (1939); Ehrenfest, Z. Phys. 45, 455 (1927); Robertson, Phys.
# Rev. 34, 163 (1929); Margolus & Levitin, Physica D 120, 188 (1998).
# (Mandelstam & Tamm, J. Phys. USSR 9, 249 (1945) predates DOIs — cited
# inline only, no bib entry.)

"""
    VirialTheorem <: AbstractRelation

The quantum virial theorem for a homogeneous potential of degree `n`
(`V(λr) = λⁿ V(r)`),

`2⟨T⟩ = n⟨V⟩`,

(Euler's theorem on the stationary state).  Harmonic `n = 2` gives
`⟨T⟩ = ⟨V⟩`; Coulomb `n = −1` gives `2⟨T⟩ = −⟨V⟩`, so `E = −⟨T⟩ = ½⟨V⟩`.

Variables: `T` = `⟨T⟩`, `V` = `⟨V⟩`, `n`.
"""
@relation :quantum VirialTheorem(T::KineticEnergy, V::PotentialEnergy, n) = 2 * T - n * V

"""
    HellmannFeynman <: AbstractRelation

The Hellmann–Feynman theorem (Feynman, Phys. Rev. 56, 340 (1939)): the
derivative of an eigenenergy with respect to a parameter is the
expectation of the derivative of the Hamiltonian,

`dE/dλ = ⟨∂H/∂λ⟩`.

Supplied-derivative convention: `dH_dλ` is the caller-computed expectation
`⟨∂H/∂λ⟩`.  Variables: `dE_dλ`, `dH_dλ`.
"""
@relation :quantum HellmannFeynman(dE_dλ, dH_dλ) = dE_dλ - dH_dλ

"""
    EhrenfestPosition <: AbstractRelation

The Ehrenfest theorem for position (Ehrenfest, Z. Phys. 45, 455 (1927)):
the mean position moves with the mean velocity,

`d⟨x⟩/dt = ⟨p⟩ / m`.

Variables: `dx_dt`, `p` = `⟨p⟩`, `m`.
"""
@relation :quantum EhrenfestPosition(dx_dt, p, m) = dx_dt - p / m

"""
    EhrenfestMomentum <: AbstractRelation

The Ehrenfest theorem for momentum (Ehrenfest, Z. Phys. 45, 455 (1927)):
the mean momentum obeys the classical force law,

`d⟨p⟩/dt = −⟨∂V/∂x⟩ = ⟨F⟩`,

the quantum counterpart of Newton's second law.  Variables: `dp_dt`,
`F` = `⟨F⟩`.
"""
@relation :quantum EhrenfestMomentum(dp_dt, F) = dp_dt - F

"""
    EnergyVarianceEigenstate <: AbstractRelation

The zero-variance eigenstate condition: an exact eigenstate has vanishing
energy variance,

`⟨H²⟩ = ⟨H⟩²`   ⟺   `Var(H) = ⟨H²⟩ − E² = 0`   (`E = ⟨H⟩`),

the standard convergence metric of a variational / DMRG ground-state
calculation (`Var(H) → 0` as the state approaches an eigenstate).

Variables: `H2` = `⟨H²⟩`, `E` = `⟨H⟩`.
"""
@relation :quantum EnergyVarianceEigenstate(H2, E) = H2 - E^2

"""
    RobertsonUncertainty <: AbstractInequality

The Robertson uncertainty relation (Robertson, Phys. Rev. 34, 163
(1929)),

`ΔA · ΔB ≥ ½ |⟨[A, B]⟩|`

(slack `ΔA·ΔB − ½|⟨[A,B]⟩|`), generalizing Heisenberg `Δx·Δp ≥ ℏ/2`
(`|⟨[x,p]⟩| = ℏ`).  Saturated by a minimum-uncertainty (coherent /
squeezed) state.

Variables: `ΔA`, `ΔB`, `comm` = `⟨[A, B]⟩`.
"""
@inequality :quantum RobertsonUncertainty(ΔA, ΔB, comm) = ΔA * ΔB - abs(comm) / 2

"""
    LiebRobinsonBound <: AbstractInequality

The Lieb–Robinson bound (Lieb & Robinson, Commun. Math. Phys. 28, 251
(1972)): information in a locally-interacting quantum system spreads no
faster than an emergent velocity `v_LR` — the group velocity of
correlations is bounded,

`v ≤ v_LR`

(slack `v_LR − v`).  An effective light cone; the many-body analogue of
relativistic causality, setting entanglement-growth and thermalization
rates.

Variables: `v`, `v_LR`.
"""
@inequality :quantum LiebRobinsonBound(v, v_LR) = v_LR - v

# ─── Quantum speed limits: the minimal time to an orthogonal state ───────

"""
    MandelstamTammBound <: AbstractInequality

The Mandelstam–Tamm quantum speed limit: the time to evolve to an orthogonal
state is bounded below by the energy uncertainty (Mandelstam & Tamm, J. Phys.
(USSR) 9, 249 (1945); `ħ = 1`),

`τ⊥ ≥ π / (2 ΔE)`,   `ΔE = √(⟨H²⟩ − ⟨H⟩²)`

(slack `τ − π/(2 ΔE)`).  The energy–time bound: no state evolves to a distinguishable
one faster than its energy spread allows.

Variables: `τ` = orthogonalization time, `ΔE` = energy uncertainty.
"""
@inequality :quantum MandelstamTammBound(τ, ΔE) = τ - π / (2 * ΔE)

"""
    MargolusLevitinBound <: AbstractInequality

The Margolus–Levitin quantum speed limit: the orthogonalization time is *also*
bounded below by the mean energy above the ground state (Margolus & Levitin,
Physica D 120, 188 (1998); `ħ = 1`),

`τ⊥ ≥ π / (2 (E − E₀))`

(slack `τ − π/(2 E_above)`, `E_above = E − E₀`).  Independent of and complementary to
[`MandelstamTammBound`](@ref); the true limit is set by whichever is tighter,
`τ⊥ ≥ (π/2) / min(ΔE, E − E₀)`.

Variables: `τ` = orthogonalization time, `E_above` = `E − E₀` (mean energy above the
ground state).
"""
@inequality :quantum MargolusLevitinBound(τ, E_above) = τ - π / (2 * E_above)
