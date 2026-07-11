# core/fields.jl — the conjugate-field vocabulary.
#
# The thermodynamic potentials are functions of a set of control fields
# (temperature, magnetic field, chemical potential, …), and every
# response function is a derivative with respect to one of them.  These
# fields are the differentiation variables of the response genealogy
# (structure/response.jl); as tags they let that genealogy — and any
# downstream code — dispatch on "with respect to what".

"""
    AbstractField

Abstract parent type for the intensive control fields a thermodynamic
potential depends on.  Concrete singleton tags — [`Temperature`](@ref),
[`InverseTemperature`](@ref), [`MagneticField`](@ref),
[`ChemicalPotential`](@ref) — name the variable a response function is
differentiated with respect to.
"""
abstract type AbstractField end
export AbstractField

"""
    Temperature <: AbstractField

The temperature `T`.  Conjugate (via `S = −∂F/∂T`) to the entropy.
"""
struct Temperature <: AbstractField end
export Temperature

"""
    InverseTemperature <: AbstractField

The inverse temperature `β = 1/T`.  The natural variable of the
Gibbs–Helmholtz relation `U = ∂(βF)/∂β`.
"""
struct InverseTemperature <: AbstractField end
export InverseTemperature

"""
    MagneticField <: AbstractField

The magnetic field `h`.  Conjugate (via `M = −∂F/∂h`) to the
magnetization — see [`conjugate_field`](@ref).
"""
struct MagneticField <: AbstractField end
export MagneticField

"""
    ChemicalPotential <: AbstractField

The chemical potential `μ`.  Conjugate to the particle number
(`N = −∂Ω/∂μ` for the grand potential Ω); the defining field of the
[`GrandCanonical`](@ref) ensemble.
"""
struct ChemicalPotential <: AbstractField end
export ChemicalPotential

"""
    conjugate_field(quantity) -> AbstractField
    conjugate_field(::Type{<:AbstractQuantity}) -> AbstractField

The field a quantity is thermodynamically conjugate to: the field whose
derivative of the free energy *is* (up to sign) that quantity.  The
magnetization is conjugate to the [`MagneticField`](@ref)
(`M = −∂F/∂h`); the entropy to the [`Temperature`](@ref).  Undefined for
quantities that are not first field-derivatives of the free energy.
"""
conjugate_field(q::AbstractQuantity) = conjugate_field(typeof(q))
conjugate_field(::Type{<:AbstractMagnetization}) = MagneticField()
conjugate_field(::Type{ThermalEntropy}) = Temperature()
export conjugate_field
