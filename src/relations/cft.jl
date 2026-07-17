# relations/cft.jl — conformal-field-theory finite-size relations.
#
# The identities that read a 1D critical point's universal data — the
# central charge and the operator scaling dimensions — straight out of
# FINITE-SIZE spectra, i.e. exactly what an MPS / ED / DMRG calculation
# produces.  All for a periodic chain of length L with velocity v.

"""
    CasimirCentralCharge <: AbstractRelation

The universal finite-size (Casimir) correction to the ground-state
energy *density* of a periodic 1D critical chain reads off the central
charge,

`e₀(L) = e_∞ − π c v / (6 L²)`

(Blöte, Cardy & Nightingale, Phys. Rev. Lett. 56, 742 (1986); Affleck,
Phys. Rev. Lett. 56, 746 (1986)).  Supplied-value convention:
`dE = e₀(L) − e_∞` is the caller-computed finite-size correction to the
ground-state energy per site.

Variables: `dE`, `c`, `v`, `L`.
"""
@relation :cft CasimirCentralCharge(dE, c::CentralCharge, v, L) = dE + π * c * v / (6 * L^2)

"""
    FiniteSizeGap <: AbstractRelation

The finite-size energy gap of a periodic 1D critical chain gives the
scaling dimension of the corresponding operator,

`E_x(L) − E₀(L) = 2π v x / L`

(Cardy, Nucl. Phys. B 270, 186 (1986)): each primary/descendant with
scaling dimension `x` appears as a level whose gap closes as `1/L` with a
universal amplitude `2πvx`.  Reads `x` off the measured finite-size gap.

Variables: `gap` = `E_x(L) − E₀(L)`, `x`, `v`, `L`.
"""
@relation :cft FiniteSizeGap(gap::MassGap, x::ScalingDimension, v, L) = gap - 2π * v * x / L

"""
    CardyDensityOfStates <: AbstractRelation

Cardy's asymptotic density of states of a 2D CFT (Cardy, Nucl. Phys. B
270, 186 (1986)): the number of states at large scaling dimension `Δ`
grows as

`ln ρ(Δ) = 2π √(c Δ / 6)`

(the modular-invariance image of the ground-state Casimir energy; fixes
the microcanonical entropy of a CFT from its central charge `c`).

Variables: `ln_ρ` = `ln ρ(Δ)`, `c`, `Δ`.
"""
@relation :cft CardyDensityOfStates(ln_ρ, c::CentralCharge, Δ::ScalingDimension) =
    ln_ρ - 2π * sqrt(c * Δ / 6)

"""
    CTheorem <: AbstractInequality

Zamolodchikov's c-theorem (Zamolodchikov, JETP Lett. 43, 730 (1986)): the
central charge decreases monotonically along renormalization-group flow
from the ultraviolet to the infrared fixed point,

`c_UV ≥ c_IR`

(slack `c_UV − c_IR`).  Irreversibility of RG flow — `c` counts the
massless degrees of freedom, which can only be integrated out.

Variables: `c_UV`, `c_IR`.
"""
@inequality :cft CTheorem(c_UV, c_IR) = c_UV - c_IR
