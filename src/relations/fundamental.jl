# relations/fundamental.jl — the fundamental equations of thermodynamics.
#
# The exact algebraic web connecting Z, F, U, S at fixed temperature:
#
#   F = −β⁻¹ ln Z           (statistical definition — micro ↔ macro bridge)
#   F = U − T·S              (Helmholtz Legendre transform)
#   S = −∂F/∂T               (entropy as a free-energy response)
#   U = ∂(βF)/∂β             (Gibbs–Helmholtz, β form)
#
# One @relation declaration each.  The derivative forms follow the
# supplied-derivative convention (cf. `LinearResponseFDT`): the caller
# provides the derivative value however obtained — closed form, AD,
# finite difference — and the relation states what it must equal.
# Units must be homogeneous: all-total or all-per-site (`N` on
# `FreeEnergyFromZ`; `F`, `U`, `S` then share that granularity).

"""
    FreeEnergyFromZ <: AbstractRelation

The statistical definition of the Helmholtz free energy,

`f = −ln(Z) / (β N)`,

bridging the microscopic partition function and the macroscopic
potential.  `N = 1` (default) gives the total free energy; `N` = number
of sites gives the per-site density (the QAtlas `FreeEnergy` tag
convention, `f = -β⁻¹ log Z / N`).  Note the log makes this relation
inherently floating-point — the exact-arithmetic contract applies only
to the arithmetic around it.  Non-affine in `Z`, so `Val(:Z)` has a
specialized solve (the exp inverse); every other variable is generic.
"""
@relation :fundamental FreeEnergyFromZ(f, Z, β, N=1) = f - (-log(Z) / (β * N))

function _solve(::FreeEnergyFromZ, ::Val{:Z}; f, β, N=1, _extra...)
    return exp(-β * N * f)
end

"""
    FreeEnergyLegendre <: AbstractRelation

The fundamental (Helmholtz–Legendre) relation among the potentials at
fixed temperature,

`F = U − T·S`   ⟺   `S = β(U − F)`,

with all three potentials in the same granularity (all-total or
all-per-site).  Purely algebraic: exact inputs give exact residuals.
"""
@relation :fundamental FreeEnergyLegendre(F, U, S, β) = F - (U - S / β)

"""
    EntropyResponse <: AbstractRelation

Entropy as a free-energy response,

`S = −∂F/∂T`.

Supplied-derivative convention: the caller provides `dF_dT` however
obtained (closed form, AD, finite difference).  Reconciling this
derivative route against the algebraic [`FreeEnergyLegendre`](@ref)
route is the classic thermodynamic self-consistency check.
"""
@relation :fundamental EntropyResponse(S, dF_dT) = S - (-dF_dT)

"""
    GibbsHelmholtz <: AbstractRelation

The Gibbs–Helmholtz equation in the β form,

`U = ∂(βF)/∂β`.

Supplied-derivative convention: `dβF_dβ` is the caller-computed value
of `∂(βF)/∂β` (equivalently `−∂ln Z/∂β`, since `βF = −ln Z`) evaluated
at the same state point as `U`.
"""
@relation :fundamental GibbsHelmholtz(U, dβF_dβ) = U - dβF_dβ

"""
    MagnetizationResponse <: AbstractRelation

The order parameter as the field-derivative of the free energy,

`M = −∂F/∂h`.

The first edge of the field-derivative genealogy
([`derivative_edge`](@ref)`(MagnetizationZ)`), stated exactly.
Supplied-derivative convention: `dF_dh` is the caller-computed
`∂F/∂h` at the working point.
"""
@relation :fundamental MagnetizationResponse(M, dF_dh) = M - (-dF_dh)

"""
    SusceptibilityResponse <: AbstractRelation

The susceptibility as the field-derivative of the order parameter,

`χ = ∂M/∂h  ( = −∂²F/∂h² )`.

The second field-derivative edge of the genealogy
([`derivative_edge`](@ref)`(SusceptibilityZZ)`), stated exactly — the
*definitional* companion of the *statistical* [`SusceptibilityFDT`](@ref)
(`χ = β·Var(M)`): the same response reached two ways.  Supplied-
derivative convention: `dM_dh` is the caller-computed `∂⟨M⟩/∂h`.
"""
@relation :fundamental SusceptibilityResponse(χ, dM_dh) = χ - dM_dh
