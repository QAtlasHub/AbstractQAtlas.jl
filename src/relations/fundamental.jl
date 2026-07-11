# relations/fundamental.jl — the fundamental equations of thermodynamics.
#
# The exact algebraic web connecting Z, F, U, S at fixed temperature:
#
#   F = −β⁻¹ ln Z           (statistical definition — micro ↔ macro bridge)
#   F = U − T·S              (Helmholtz Legendre transform)
#   S = −∂F/∂T               (entropy as a free-energy response)
#   U = ∂(βF)/∂β             (Gibbs–Helmholtz, β form)
#
# The first two are purely algebraic and close under the three-verb
# interface directly.  The derivative forms follow the same convention
# as `LinearResponseFDT`: the caller supplies the derivative *value*
# (however obtained — closed form, AD, finite difference), and the
# relation states what it must equal.  Comparing an independently
# computed derivative route against the algebraic route is exactly the
# reconciliation the QAtlas identities plane performs.
#
# Units must be homogeneous: either all-total or all-per-site (divide
# ln Z by N via the `N` variable of `FreeEnergyFromZ`; `F`, `U`, `S`
# then all carry the same granularity).

"""
    FreeEnergyFromZ <: AbstractRelation

The statistical definition of the Helmholtz free energy,

`f = −ln(Z) / (β N)`,

bridging the microscopic partition function and the macroscopic
potential.  `N = 1` (default) gives the total free energy; `N` = number
of sites gives the per-site density (the QAtlas `FreeEnergy` tag
convention, `f = -β⁻¹ log Z / N`).

Variables: `f`, `Z`, `β` (or `T`), `N` (default 1).  Note the log makes
this relation inherently floating-point — the exact-arithmetic contract
applies only to the arithmetic around it.
"""
struct FreeEnergyFromZ <: AbstractRelation end
export FreeEnergyFromZ

function residual(::FreeEnergyFromZ; f, Z, β=nothing, T=nothing, N=1)
    b = _beta(; β=β, T=T)
    return f - (-log(Z) / (b * N))
end
function solve(::FreeEnergyFromZ, ::Val{:f}; Z, β=nothing, T=nothing, N=1)
    b = _beta(; β=β, T=T)
    return -log(Z) / (b * N)
end
function solve(::FreeEnergyFromZ, ::Val{:Z}; f, β=nothing, T=nothing, N=1)
    b = _beta(; β=β, T=T)
    return exp(-b * N * f)
end

"""
    FreeEnergyLegendre <: AbstractRelation

The fundamental (Helmholtz–Legendre) relation among the potentials at
fixed temperature,

`F = U − T·S`   ⟺   `S = β(U − F)`,

with all three potentials in the same granularity (all-total or
all-per-site).  Purely algebraic: exact inputs give exact residuals.

Variables: `F`, `U`, `S`, `β` (or `T`).
"""
struct FreeEnergyLegendre <: AbstractRelation end
export FreeEnergyLegendre

function residual(::FreeEnergyLegendre; F, U, S, β=nothing, T=nothing)
    b = _beta(; β=β, T=T)
    return F - (U - S / b)
end
function solve(::FreeEnergyLegendre, ::Val{:F}; U, S, β=nothing, T=nothing)
    b = _beta(; β=β, T=T)
    return U - S / b
end
function solve(::FreeEnergyLegendre, ::Val{:U}; F, S, β=nothing, T=nothing)
    b = _beta(; β=β, T=T)
    return F + S / b
end
function solve(::FreeEnergyLegendre, ::Val{:S}; F, U, β=nothing, T=nothing)
    b = _beta(; β=β, T=T)
    return b * (U - F)
end

"""
    EntropyResponse <: AbstractRelation

Entropy as a free-energy response,

`S = −∂F/∂T`.

Following the supplied-derivative convention (cf.
[`LinearResponseFDT`](@ref)): the caller provides the derivative value
`dF_dT` however obtained (closed form, AD, finite difference), and the
relation states what it must equal.  Reconciling this derivative route
against the algebraic [`FreeEnergyLegendre`](@ref) route is the classic
thermodynamic self-consistency check.

Variables: `S`, `dF_dT`.
"""
struct EntropyResponse <: AbstractRelation end
export EntropyResponse

residual(::EntropyResponse; S, dF_dT) = S - (-dF_dT)
solve(::EntropyResponse, ::Val{:S}; dF_dT) = -dF_dT
solve(::EntropyResponse, ::Val{:dF_dT}; S) = -S

"""
    GibbsHelmholtz <: AbstractRelation

The Gibbs–Helmholtz equation in the β form,

`U = ∂(βF)/∂β`.

Supplied-derivative convention: `dβF_dβ` is the caller-computed value
of `∂(βF)/∂β` (equivalently `−∂ln Z/∂β`, since `βF = −ln Z`) evaluated
at the same state point as `U`.

Variables: `U`, `dβF_dβ`.
"""
struct GibbsHelmholtz <: AbstractRelation end
export GibbsHelmholtz

residual(::GibbsHelmholtz; U, dβF_dβ) = U - dβF_dβ
solve(::GibbsHelmholtz, ::Val{:U}; dβF_dβ) = dβF_dβ
solve(::GibbsHelmholtz, ::Val{:dβF_dβ}; U) = U
