# relations/ensembles.jl — statistical-ensemble relations and the
# thermal-pure-quantum-state estimators.
#
# The bridge between the microcanonical and canonical descriptions, plus
# the thermal-pure-quantum (TPQ) formulas that realize them from a single
# random vector — the kind of identity a many-body calculation (ED, MPS,
# TPQ) checks its measured quantities against.
#
# Ensemble web:
#   microcanonical  S(E),  β = ∂S/∂E            (MicrocanonicalTemperature)
#        ⟷ canonical via the Legendre transform  S(E) = β(E − F)
#          (that Legendre step is `FreeEnergyLegendre` in fundamental.jl),
#          the two ensembles agreeing at E = U(β) in the thermodynamic
#          limit (ensemble equivalence).
#   canonical Z(β) from a random state          (CanonicalTPQ)
#     Sugiura & Shimizu, Phys. Rev. Lett. 108, 240401 (2012) [micro];
#                        Phys. Rev. Lett. 111, 010401 (2013) [canonical].

"""
    MicrocanonicalTemperature <: AbstractRelation

The microcanonical (Boltzmann–Gibbs) definition of the inverse
temperature as the energy-derivative of the entropy,

`β = ∂S/∂E`

with `S(E)` the microcanonical entropy (`S = ln W(E)`, `W` the number of
states in the energy shell).  Supplied-derivative convention: `dS_dE` is
the caller-computed `∂S/∂E` at the working energy.  Ensemble equivalence
identifies this microcanonical `β` with the canonical control parameter
at `E = U(β)` — the connection a finite-`T` calculation can cross-check.

Variables: `β`, `dS_dE`.
"""
@relation :ensemble MicrocanonicalTemperature(β::InverseTemperature, dS_dE) = β - dS_dE

"""
    CanonicalTPQ <: AbstractRelation

The canonical thermal-pure-quantum estimator of the partition function
(Sugiura & Shimizu, Phys. Rev. Lett. 111, 010401 (2013)),

`Z(β) = D · ⟨ψ₀| e^{−βĤ} |ψ₀⟩`,

where `|ψ₀⟩` is a Haar-random normalized state in the `D`-dimensional
Hilbert space and the bar is the random-state average — exact because
`⟨ψ₀| Ô |ψ₀⟩ = Tr Ô / D` on average, with fluctuations exponentially
small in system size.  Thermal averages follow the same way,
`⟨Â⟩_β = ⟨ψ_β| Â |ψ_β⟩ / ⟨ψ_β|ψ_β⟩` with `|ψ_β⟩ = e^{−βĤ/2}|ψ₀⟩`.

Supplied-weight convention: `tpq_weight = ⟨ψ₀| e^{−βĤ} |ψ₀⟩`.
Variables: `Z`, `tpq_weight`, `D`.
"""
@relation :ensemble CanonicalTPQ(Z::PartitionFunction, tpq_weight, D) = Z - D * tpq_weight
