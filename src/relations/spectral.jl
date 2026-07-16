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
@relation :spectral Dyson(G, G0, Σ) = inv(G) - (inv(G0) - Σ)

"""
    SpectralFromGreens <: AbstractRelation

The spectral representation of the retarded Green's function at `(q, ω)`:

`A = −(1/π) Im G^R`   ⟺   `A + Im(G^R)/π = 0`.

Takes the **complex** retarded Green's function `G = G^R(q, ω)` directly
(the same variable Dyson speaks of, so one `G` in a measurement bag fires
both), and the real spectral weight `A`; only `Im G^R` enters here (the
real part is fixed separately by Kramers–Kronig).  This makes the relation
turnkey against a raw measurement — no need to pre-project the imaginary
part.  Exact-arithmetic-safe: `imag` of a `Complex{Rational}` stays
`Rational`.
"""
@relation :spectral SpectralFromGreens(A, G) = A + imag(G) / π

"""
    SpectralSumRule <: AbstractRelation

The single-band spectral normalization

`∫ A(q, ω) dω = 1`.

Supplied-integral convention: `spectral_integral` is the caller-computed
frequency integral of the spectral function at fixed `q`.
"""
@relation :spectral SpectralSumRule(spectral_integral) = spectral_integral - 1

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
@relation :spectral NMRExponent(θ_NMR, Δ_op) = θ_NMR - (2 * Δ_op - 1)

"""
    StaticFromDynamicalStructureFactor <: AbstractRelation

The static (equal-time) structure factor as the frequency integral of
the dynamical one,

`S(q) = ∫ S(q, ω) dω / (2π)`

(Van Hove, Phys. Rev. 95, 249 (1954)).  Supplied-integral convention:
`sqw_integral = ∫ S(q, ω) dω/(2π)` is the caller-computed frequency
integral at fixed `q`.

Variables: `Sq`, `sqw_integral`.
"""
@relation :spectral StaticFromDynamicalStructureFactor(Sq, sqw_integral) = Sq - sqw_integral

"""
    DynamicalFDT <: AbstractRelation

The finite-temperature fluctuation–dissipation theorem relating the
dynamical structure factor to the dissipative part of the dynamical
susceptibility,

`S(q, ω) = χ''(q, ω) / [π (1 − e^{−βω})]`

(Callen & Welton, Phys. Rev. 83, 34 (1951)).  Since `χ''` is odd in `ω`,
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
@relation :spectral CorrelationLengthGap(ξ, v, Δ) = ξ - v / Δ

"""
    KramersKronigReal <: AbstractRelation

The Kramers–Kronig relation fixing the **real** part of a causal
(retarded, analytic-in-the-upper-half-plane) response function from the
Hilbert transform of its imaginary part,

`χ'(ω) = (1/π) P ∫ χ''(ω') / (ω' − ω) dω'`,

(Kronig, J. Opt. Soc. Am. 12, 547 (1926); Toll, Phys. Rev. 104, 1760
(1956)).  Applies to any causal response — the optical conductivity
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

the companion of [`KramersKronigReal`](@ref) (Kronig, J. Opt. Soc. Am. 12,
547 (1926); Toll, Phys. Rev. 104, 1760 (1956)).  Supplied-integral
convention: `pv_real` is the caller-computed principal-value Hilbert
transform `P ∫ χ'(ω')/(ω' − ω) dω'`.

Variables: `Imχ` = `χ''(ω)`, `pv_real`.
"""
@relation :spectral KramersKronigImag(Imχ, pv_real) = Imχ + pv_real / π
