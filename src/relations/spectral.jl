# relations/spectral.jl — the pointwise identities relating the
# dynamical / spectral quantities.
#
# These are the frequency-resolved analogues of the thermodynamic
# identities: each holds at a single (q, ω) point (or, for the sum rule,
# given a supplied frequency integral), so each is a first-class
# `@relation` checkable against numbers.  The transform / BZ-sum / low-ω-
# limit relations that are NOT single-point live in `structure/spectral.jl`
# as the inter-quantity genealogy.
#
# Domain tag :spectral throughout, so `relation_report` can sweep a bag of
# measured propagator/response values and report every identity that
# applies.

"""
    Dyson <: AbstractRelation

The Dyson equation relating the full and bare propagators through the
self-energy, at fixed `(q, ω)`:

`G^{-1} = G₀^{-1} − Σ`.

Written with `inv` rather than `1/…`, so the SAME identity is honest
about the orbital-tensor character of the propagators: it holds
verbatim for scalar (single-band) `G, G0, Σ` AND for matrix-valued
`G_ab, Σ_ab` in orbital/band space — `residual` returns the residual
matrix, whose norm should vanish.  (`check`/`solve` are scalar; matrix
inputs use `residual` + a norm.)  Complex-valued.
"""
@relation :spectral Dyson(G::RetardedGreensFunction, G0, Σ::SelfEnergy) =
    inv(G) - (inv(G0) - Σ)

"""
    SpectralFromGreens <: AbstractRelation

The spectral representation of the retarded Green's function at `(q, ω)`:

`A = −(1/π) Im G^R`   ⟺   `A + Im(G^R)/π = 0`.

Pass the full (complex) retarded Green's function `G = G^R(q, ω)` and the real
spectral weight `A`; the imaginary part is taken in the kernel.  Keyed on the
`RetardedGreensFunction` TYPE — the same `G` as [`Dyson`](@ref) and the Keldysh
relations, never a separate `ImGR`/`GR` symbol.
"""
@relation :spectral SpectralFromGreens(A::SpectralFunction, G::RetardedGreensFunction) =
    A + imag(G) / π

"""
    SpectralSumRule <: AbstractRelation

The single-band spectral normalization

`∫ A(q, ω) dω = 1`.

Supplied-integral convention: `spectral_integral` is the caller-computed
frequency integral of the spectral function at fixed `q`.
"""
@relation :spectral SpectralSumRule(spectral_integral) = spectral_integral - 1

"""
    FSumRule <: AbstractRelation

The **f-sum rule** — the first frequency moment of the dynamical structure factor,

`∫ ω S(q, ω) dω = N q² / (2m)`

(`ℏ = 1`; `N` particles of mass `m`, with `N = 1` the per-particle form).  A
model-independent identity: the first moment equals `½⟨[[H, ρ_q], ρ_{−q}]⟩`, and for a
`q²/2m` kinetic energy that double commutator is `N q²/m` regardless of the
interactions (Pines & Nozières, *The Theory of Quantum Liquids* 1966; Thomas–Reiche–Kuhn
1925).  Supplied-integral convention: `first_moment = ∫ ω S(q, ω) dω` is the
caller-computed first moment at fixed `q`.

Variables: `first_moment`, `q`, `m`, `N = 1`.
"""
@relation :spectral FSumRule(first_moment, q, m, N=1) = first_moment - N * q^2 / (2 * m)

"""
    DetailedBalance <: AbstractRelation

The finite-temperature detailed-balance condition on the dynamical
structure factor,

`S(q, −ω) = e^{−βω} S(q, ω)`,

a convention-independent consequence of the fluctuation–dissipation
theorem.  Variables: `S_plus = S(q, ω)`, `S_minus = S(q, −ω)`, `ω`, and
`β` (or `T`).
"""
@relation :spectral DetailedBalance(S_plus, S_minus, ω, β) = S_minus - exp(-β * ω) * S_plus

"""
    NMRExponent <: AbstractRelation

The dynamical scaling relation fixing the NMR spin–lattice relaxation
exponent from the operator scaling dimension at a quantum critical
point,

`θ_NMR = 2 Δ_op − 1`   (with `1/T₁ ∝ T^{θ_NMR}` as `T → 0`).

Exact-arithmetic: `Δ_op = 1//8` (1D TFIM QCP) gives `θ_NMR = −3//4`
exactly.
"""
@relation :spectral NMRExponent(θ_NMR::NMRRelaxationExponent, Δ_op::ScalingDimension) =
    θ_NMR - (2 * Δ_op - 1)

"""
    StaticFromDynamicalStructureFactor <: AbstractRelation

The static (equal-time) structure factor as the frequency integral of
the dynamical one,

`S(q) = ∫ S(q, ω) dω / (2π)`

(Van Hove, [VanHove1954](@cite)).  Supplied-integral convention:
`sqw_integral = ∫ S(q, ω) dω/(2π)` is the caller-computed frequency
integral at fixed `q`.

Variables: `Sq`, `sqw_integral`.
"""
@relation :spectral StaticFromDynamicalStructureFactor(
    Sq::StaticStructureFactor, sqw_integral
) = Sq - sqw_integral

"""
    StaticStructureFactorFromCorrelation <: AbstractRelation

The static structure factor at zero wavevector as the spatial integral of the
two-point correlation — the static fluctuation / compressibility sum rule,

`S(q → 0) = ∫ G(r) dr`

(Chaikin–Lubensky, *Principles of Condensed Matter Physics*).  Supplied-integral
convention: `integral_G = ∫ G(r) dr` is the caller-computed spatial integral of the
correlation function (evaluating it is the functional sibling's job, issue #14).
Together with [`StructureFactorSusceptibility`](@ref) (`χ = β S(q→0)`) this closes the
static loop `χ = ∫G(r)dr = S(q→0)`; the structural edge is
`spectral_origin(StaticStructureFactor)`.

Variables: `Sq0`, `integral_G`.
"""
@relation :spectral StaticStructureFactorFromCorrelation(
    Sq0::StaticStructureFactor, integral_G
) = Sq0 - integral_G

"""
    DynamicalFDT <: AbstractRelation

The finite-temperature fluctuation–dissipation theorem relating the
dynamical structure factor to the dissipative part of the dynamical
susceptibility,

`S(q, ω) = χ''(q, ω) / [π (1 − e^{−βω})]`

(Callen & Welton, [CallenWelton1951](@cite)).  Since `χ''` is odd in `ω`,
this convention reproduces detailed balance `S(q,−ω) = e^{−βω} S(q, ω)`
([`DetailedBalance`](@ref)) automatically.

Variables: `S` = `S(q, ω)`, `χpp` = `χ''(q, ω)`, `ω`, and `β` (or `T`).
"""
@relation :spectral DynamicalFDT(S, χpp, ω, β) = S - χpp / (π * (1 - exp(-β * ω)))

"""
    CorrelationLengthGap <: AbstractRelation

The correlation length of a gapped phase set by the gap and velocity,

`ξ = v / Δ`,

the real-space decay length of a relativistic dispersion
`E(k) = √(Δ² + v²k²)` (`⟨O(r)O(0)⟩ ∼ e^{−r/ξ}`).  A staple consistency
check for a gapped MPS/DMRG calculation: the measured correlation length
and the measured gap must satisfy `ξΔ = v`.

Variables: `ξ`, `v`, `Δ`.
"""
@relation :spectral CorrelationLengthGap(ξ::CorrelationLength, v::Velocity, Δ::MassGap) =
    ξ - v / Δ

"""
    MassGapPositivity <: AbstractInequality

The spectral gap is non-negative,

`Δ ≥ 0`

(slack `Δ`), because `Δ = E₁ − E₀` and `E₀` is by definition the lowest level —
so this is airtight rather than a modelling assumption.  Saturated by any
gapless phase.  A measured `Δ < 0` means the reported "ground state" was not the
ground state: the variational state is above a level the solver missed, or two
levels were ordered wrongly.

Stated on [`MassGap`](@ref) alone, NOT on the [`AbstractGap`](@ref) group.  The
sector-resolved gaps do not share it: [`SpinGap`](@ref) `= E₀(Sᶻ=1) − E₀(Sᶻ=0)`
is negative whenever the ground state is polarised (a ferromagnet), and
[`ChargeGap`](@ref) `= E₀(N+1) + E₀(N−1) − 2E₀(N)` is negative where `E₀(N)` is
concave (phase separation).  Both are excitation energies only when the
reference sector holds the ground state, which is a statement about the model.

Variables: `Δ`.
"""
@bound :spectral MassGapPositivity(Δ::MassGap >= 0)

"""
    KramersKronigReal <: AbstractRelation

The Kramers–Kronig relation fixing the **real** part of a causal
(retarded, analytic-in-the-upper-half-plane) response function from the
Hilbert transform of its imaginary part,

`χ'(ω) = (1/π) P ∫ χ''(ω') / (ω' − ω) dω'`,

(Kronig, [Kronig1926](@cite); Toll, [Toll1956](@cite)).  Applies to any causal response — the optical conductivity
`σ(ω)`, the susceptibility `χ(ω)`, the retarded Green's function, the
dielectric function `ε(ω)`.  Supplied-integral convention: `pv_imag` is
the caller-computed principal-value Hilbert transform
`P ∫ χ''(ω')/(ω' − ω) dω'`.

Variables: `Reχ` = `χ'(ω)`, `pv_imag`.
"""
@relation :spectral KramersKronigReal(Reχ, pv_imag) = Reχ - pv_imag / π

"""
    KramersKronigImag <: AbstractRelation

The Kramers–Kronig relation fixing the **imaginary** part of a causal
response function from the Hilbert transform of its real part,

`χ''(ω) = −(1/π) P ∫ χ'(ω') / (ω' − ω) dω'`,

the companion of [`KramersKronigReal`](@ref) (Kronig, [Kronig1926](@cite); Toll, [Toll1956](@cite)).  Supplied-integral
convention: `pv_real` is the caller-computed principal-value Hilbert
transform `P ∫ χ'(ω')/(ω' − ω) dω'`.

Variables: `Imχ` = `χ''(ω)`, `pv_real`.
"""
@relation :spectral KramersKronigImag(Imχ, pv_real) = Imχ + pv_real / π

"""
    ResponseRealityReal <: AbstractRelation

The **reality** (parity) condition on the real part of a causal response.  The response
of a real observable to a real field is real in time, so its Fourier transform is
conjugate-symmetric under negating all frequencies, `χ⁽ⁿ⁾(−ω⃗) = χ⁽ⁿ⁾(ω⃗)*` — hence
**`Re χ` is EVEN**:

`Re χ⁽ⁿ⁾(−ω⃗) = Re χ⁽ⁿ⁾(ω⃗)`.

It holds at every order with the same shape (linear: `Re χ(−ω) = Re χ(ω)`).  With
[`ResponseRealityImag`](@ref) (`Im χ` odd) and [`intrinsic_permutation_symmetric`](@ref)
(frequency exchange) this closes the model-independent symmetry web of the multi-time
response functions (Kubo, J. Phys. Soc. Jpn. 12, 570 (1957)).  Pass
`Re_plus = Re χ⁽ⁿ⁾(ω⃗)` and `Re_minus = Re χ⁽ⁿ⁾(−ω⃗)`.

Variables: `Re_plus`, `Re_minus`.
"""
@relation :spectral ResponseRealityReal(Re_plus, Re_minus) = Re_minus - Re_plus

"""
    ResponseRealityImag <: AbstractRelation

The **reality** (parity) condition on the imaginary part of a causal response — the
companion of [`ResponseRealityReal`](@ref).  From `χ⁽ⁿ⁾(−ω⃗) = χ⁽ⁿ⁾(ω⃗)*`, **`Im χ` is
ODD**:

`Im χ⁽ⁿ⁾(−ω⃗) = −Im χ⁽ⁿ⁾(ω⃗)`

(linear: `Im χ(−ω) = −Im χ(ω)` — the dissipative part, the oddness the
[`DetailedBalance`](@ref)/FDT convention already relies on).  Pass
`Im_plus = Im χ⁽ⁿ⁾(ω⃗)` and `Im_minus = Im χ⁽ⁿ⁾(−ω⃗)`.

Variables: `Im_plus`, `Im_minus`.
"""
@relation :spectral ResponseRealityImag(Im_plus, Im_minus) = Im_minus + Im_plus
