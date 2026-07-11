# core/quantities.jl — the model-independent quantity vocabulary.
#
# Ported from QAtlas `src/core/quantities.jl`: the abstract quantity
# hierarchy plus a CURATED set of stable concrete tags (the ones the
# relations layer and its first consumers reference).  The full tag set
# (~60 structs: NMR, Loschmidt/DQPT, entanglement family, structure
# factors, correlation families, velocities, CFT data, …) migrates here
# incrementally as QAtlas adopts this package — tracked upstream; do not
# bulk-copy it ahead of need.
#
# Three tags that were previously stranded inside a *model* file
# (`PartitionFunction`, `CriticalTemperature`, `SpontaneousMagnetization`
# in QAtlas' IsingSquare.jl) are given their generic home here — they are
# meaningful for any statistical-mechanics model.

# ─── Abstract quantity hierarchy ────────────────────────────────────────

"""
    AbstractThermalPotential <: AbstractQuantity

Scalar thermodynamic potentials and their densities (energy, free
energy, entropy, specific heat, …).
"""
abstract type AbstractThermalPotential <: AbstractQuantity end

"""
    AbstractMagnetization <: AbstractQuantity

Order-parameter expectation values ⟨M_α⟩ and their site-resolved
variants.
"""
abstract type AbstractMagnetization <: AbstractQuantity end

"""
    AbstractSusceptibility <: AbstractQuantity

Static linear-response susceptibilities χ_αα.
"""
abstract type AbstractSusceptibility <: AbstractQuantity end

"""
    AbstractTwoPointCorrelation <: AbstractQuantity

Two-point correlation functions ⟨O_i O_j⟩ (connected or full).
"""
abstract type AbstractTwoPointCorrelation <: AbstractQuantity end

"""
    AbstractStructureFactor <: AbstractQuantity

Momentum-space structure factors S(q).
"""
abstract type AbstractStructureFactor <: AbstractQuantity end

"""
    AbstractGap <: AbstractQuantity

Spectral gaps (mass, charge, spin channels).
"""
abstract type AbstractGap <: AbstractQuantity end

"""
    AbstractVelocity <: AbstractQuantity

Characteristic velocities (Fermi, Luttinger, Lieb-Robinson, …).
"""
abstract type AbstractVelocity <: AbstractQuantity end

"""
    AbstractEntanglementMeasure <: AbstractQuantity

Entanglement measures (von Neumann / Rényi entropies, negativity,
mutual information, …).
"""
abstract type AbstractEntanglementMeasure <: AbstractQuantity end

export AbstractThermalPotential, AbstractMagnetization, AbstractSusceptibility
export AbstractTwoPointCorrelation, AbstractStructureFactor, AbstractGap
export AbstractVelocity, AbstractEntanglementMeasure

# ─── component trait ────────────────────────────────────────────────────

"""
    component(q) -> Union{Symbol,Nothing}
    component(::Type{<:AbstractQuantity}) -> Union{Symbol,Nothing}

The component / index that a family leaf's *type name* encodes: the spin
axis of a magnetization (`:x`/`:y`/`:z`), the diagonal axis pair of a
susceptibility (`:xx`/`:yy`/`:zz`), or the excitation channel of a gap.
`nothing` for quantities that carry no component (the default).

Identities that hold per-component (e.g. the static FDT
`χ_αα = β·Var(M_α)/N`, or SU(2) isotropy `χ_xx = χ_yy = χ_zz`) pair
family members by matching `component`.
"""
component(q::AbstractQuantity) = component(typeof(q))
component(::Type{<:AbstractQuantity}) = nothing
export component

# ─── Scalar thermodynamics ──────────────────────────────────────────────

"""
    Energy{G}() <: AbstractThermalPotential
    Energy()                 # G = :natural — model-and-BC-natural granularity
    Energy(:total)           # explicit ⟨H⟩
    Energy(:per_site)        # explicit ⟨H⟩ / N

Ground-state / thermal energy expectation.  The type parameter `G` makes
the granularity (total vs per-site) a dispatch axis instead of a hidden
docstring contract.

`Energy()` resolves to the model's native granularity via the
[`native_energy_granularity`](@ref) trait.  Use the explicit
constructors when the caller needs a specific granularity.
"""
struct Energy{G} <: AbstractThermalPotential
    function Energy{G}() where {G}
        G isa Symbol || error("Energy granularity must be a Symbol, got $(typeof(G))")
        G in (:natural, :total, :per_site) ||
            error("unknown Energy granularity :$G; expected :natural, :total, or :per_site")
        return new{G}()
    end
end
Energy() = Energy{:natural}()
Energy(g::Symbol) = Energy{g}()
export Energy

"""
    native_energy_granularity(model, bc) -> :total | :per_site

Trait declaring which granularity the given `model` returns natively for
[`Energy`](@ref) at boundary condition `bc`.  Every model that supports
`Energy` must add a method per supported BC.  A missing method is caught
at the call site as a `MethodError`, which is intentional: it forces new
models to declare the convention rather than silently inheriting an
unrelated default.
"""
function native_energy_granularity end
export native_energy_granularity

"""
    FreeEnergy() <: AbstractThermalPotential

Helmholtz free energy per site, `f = -β⁻¹ log Z / N`.
"""
struct FreeEnergy <: AbstractThermalPotential end
export FreeEnergy

"""
    SpecificHeat() <: AbstractThermalPotential

Specific heat per site, `c_v(β) = β² (⟨H²⟩ − ⟨H⟩²) / N`.

The defining identity is available as a first-class relation:
[`SpecificHeatFDT`](@ref).
"""
struct SpecificHeat <: AbstractThermalPotential end
export SpecificHeat

"""
    ThermalEntropy() <: AbstractThermalPotential

Thermal entropy per site, `s = β(ε − f)`.
"""
struct ThermalEntropy <: AbstractThermalPotential end
export ThermalEntropy

"""
    PartitionFunction() <: AbstractThermalPotential

The partition function `Z(β) = Σ exp(-βE)` itself (finite systems).

Generic home for a tag that previously lived inside a model file:
any statistical-mechanics model with a finite configuration space can
register it.
"""
struct PartitionFunction <: AbstractThermalPotential end
export PartitionFunction

# ─── Order parameters and responses ─────────────────────────────────────

"""
    MagnetizationX() <: AbstractMagnetization

Uniform magnetization ⟨M_x⟩ per site.
"""
struct MagnetizationX <: AbstractMagnetization end
component(::Type{MagnetizationX}) = :x
export MagnetizationX

"""
    MagnetizationY() <: AbstractMagnetization

Uniform magnetization ⟨M_y⟩ per site.
"""
struct MagnetizationY <: AbstractMagnetization end
component(::Type{MagnetizationY}) = :y
export MagnetizationY

"""
    MagnetizationZ() <: AbstractMagnetization

Uniform magnetization ⟨M_z⟩ per site.
"""
struct MagnetizationZ <: AbstractMagnetization end
component(::Type{MagnetizationZ}) = :z
export MagnetizationZ

"""
    SpontaneousMagnetization() <: AbstractMagnetization

Spontaneous (symmetry-broken) order parameter `M(T)` in the ordered
phase; identically zero above the critical temperature.

Generic home for a tag that previously lived inside a model file.
"""
struct SpontaneousMagnetization <: AbstractMagnetization end
export SpontaneousMagnetization

"""
    SusceptibilityXX() <: AbstractSusceptibility

Static susceptibility χ_xx per site.  Defining identity:
[`SusceptibilityFDT`](@ref).
"""
struct SusceptibilityXX <: AbstractSusceptibility end
component(::Type{SusceptibilityXX}) = :xx
export SusceptibilityXX

"""
    SusceptibilityYY() <: AbstractSusceptibility

Static susceptibility χ_yy per site.
"""
struct SusceptibilityYY <: AbstractSusceptibility end
component(::Type{SusceptibilityYY}) = :yy
export SusceptibilityYY

"""
    SusceptibilityZZ() <: AbstractSusceptibility

Static susceptibility χ_zz per site.
"""
struct SusceptibilityZZ <: AbstractSusceptibility end
component(::Type{SusceptibilityZZ}) = :zz
export SusceptibilityZZ

# ─── Two-point correlations ─────────────────────────────────────────────

"""
    XXCorrelation() <: AbstractTwoPointCorrelation

Two-point correlation `⟨S^x_i S^x_j⟩` (connected or full per the
implementing atlas).  At criticality its decay is governed by the
anomalous dimension η — see the `correlation_decay` correspondence.
"""
struct XXCorrelation <: AbstractTwoPointCorrelation end
component(::Type{XXCorrelation}) = :xx
export XXCorrelation

"""
    YYCorrelation() <: AbstractTwoPointCorrelation

Two-point correlation `⟨S^y_i S^y_j⟩`.
"""
struct YYCorrelation <: AbstractTwoPointCorrelation end
component(::Type{YYCorrelation}) = :yy
export YYCorrelation

"""
    ZZCorrelation() <: AbstractTwoPointCorrelation

Two-point correlation `⟨S^z_i S^z_j⟩`.
"""
struct ZZCorrelation <: AbstractTwoPointCorrelation end
component(::Type{ZZCorrelation}) = :zz
export ZZCorrelation

# ─── Criticality ────────────────────────────────────────────────────────

"""
    CriticalTemperature() <: AbstractQuantity

Critical temperature `T_c` of a finite-temperature phase transition.

Generic home for a tag that previously lived inside a model file.
"""
struct CriticalTemperature <: AbstractQuantity end
export CriticalTemperature

"""
    CorrelationLength() <: AbstractQuantity

Correlation length ξ (units of lattice spacing).
"""
struct CorrelationLength <: AbstractQuantity end
export CorrelationLength

"""
    UniversalityClass() <: AbstractQuantity

The universality class a model's transition belongs to (returned as a
[`Universality`](@ref) tag by implementing atlases).
"""
struct UniversalityClass <: AbstractQuantity end
export UniversalityClass

"""
    CentralCharge() <: AbstractQuantity

Central charge `c` of the critical theory's CFT.
"""
struct CentralCharge <: AbstractQuantity end
export CentralCharge

# ─── Topology ───────────────────────────────────────────────────────────

"""
    TopologicalInvariant() <: AbstractQuantity

The model's topological invariant (winding number, Chern number, ℤ₂
index, Pfaffian sign, … — the concrete meaning is declared by the
implementing model).  Generic *computations* of standard invariants on
Bloch maps live in this package's relations layer
([`winding_number`](@ref), [`chern_number`](@ref)).
"""
struct TopologicalInvariant <: AbstractQuantity end
export TopologicalInvariant
