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

Linear and nonlinear response susceptibilities `χ⁽ⁿ⁾_{α;β₁…βₙ}` — see
[`Susceptibility`](@ref) for the arbitrary-order tensor.
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

# The tensor traits (`tensor_rank`, `index_spaces`, `indices`) — the
# honest successors of the old fused `component` label — live in
# `core/indices.jl`; concrete tensor quantities add their methods below.

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

Specific heat (per site) at constant volume, `c_v(β) = β² (⟨H²⟩ − ⟨H⟩²) / N`.

Defining identities: the fluctuation form [`SpecificHeatFDT`](@ref)
`c_v = β² Var(E) / N`, the entropy form [`SpecificHeatFromEntropy`](@ref)
`c_v = T ∂s/∂T`, and the difference from the constant-pressure heat
capacity [`HeatCapacityDifference`](@ref).
"""
struct SpecificHeat <: AbstractThermalPotential end
export SpecificHeat

"""
    IsobaricSpecificHeat() <: AbstractThermalPotential

Specific heat at constant pressure, `c_p`.  Exceeds the constant-volume
[`SpecificHeat`](@ref) `c_v` by `c_p − c_v = T v α² / κ_T`
([`HeatCapacityDifference`](@ref)).
"""
struct IsobaricSpecificHeat <: AbstractThermalPotential end
export IsobaricSpecificHeat

"""
    ThermalExpansionCoefficient() <: AbstractQuantity

The (isobaric, volumetric) thermal-expansion coefficient
`α = (1/V)(∂V/∂T)_p`.
"""
struct ThermalExpansionCoefficient <: AbstractQuantity end
export ThermalExpansionCoefficient

"""
    IsothermalCompressibility() <: AbstractQuantity

The isothermal compressibility `κ_T = −(1/V)(∂V/∂p)_T`.
"""
struct IsothermalCompressibility <: AbstractQuantity end
export IsothermalCompressibility

"""
    Pressure() <: AbstractQuantity

The pressure `p = −(∂F/∂V)_T` — conjugate to the [`Volume`](@ref).
"""
struct Pressure <: AbstractQuantity end
export Pressure

"""
    Volume() <: AbstractQuantity

The volume `V` — conjugate to the [`Pressure`](@ref).
"""
struct Volume <: AbstractQuantity end
export Volume

"""
    ParticleNumber() <: AbstractQuantity

The particle number `N = −(∂Ω/∂μ)_{T,V}` — conjugate to the chemical
potential ([`ChemicalPotential`](@ref)).
"""
struct ParticleNumber <: AbstractQuantity end
export ParticleNumber

"""
    LatentHeat() <: AbstractQuantity

The latent heat `L = T ΔS` of a first-order transition — the entropy
jump across the phase boundary times the temperature.  Enters the
Clausius–Clapeyron relation ([`ClausiusClapeyron`](@ref)).
"""
struct LatentHeat <: AbstractThermalPotential end
export LatentHeat

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

# ─── Order parameters and responses (tensors in spin space) ─────────────
#
# These were the quantities blurred into scalars by baking a component
# into the type name (`MagnetizationX`, `SusceptibilityZZ`).  The honest
# form carries the spin-axis index/indices as type parameters, so a
# component is a *selection* of the tensor and off-diagonal components
# (χ_xy) are expressible — not only the diagonal.  `_axis` validates that
# a parameter is a bare spin-axis symbol.

function _axis(a)
    return if a isa Symbol
        a
    else
        error("spin-axis index must be a Symbol (e.g. :x), got $(repr(a))")
    end
end

"""
    Magnetization{A}() <: AbstractMagnetization
    Magnetization(a::Symbol)

Uniform magnetization component `⟨M_A⟩` per site — a rank-1 tensor in
[`SpinAxis`](@ref) space, `A ∈ {:x, :y, :z, …}`.  `Magnetization(:z)`
replaces the old `MagnetizationZ`.
"""
struct Magnetization{A} <: AbstractMagnetization
    Magnetization{A}() where {A} = (_axis(A); new{A}())
end
Magnetization(a::Symbol) = Magnetization{a}()
tensor_rank(::Type{<:Magnetization}) = 1
index_spaces(::Type{<:Magnetization}) = (SpinAxis(),)
indices(::Type{Magnetization{A}}) where {A} = (A,)
export Magnetization

"""
    SpontaneousMagnetization() <: AbstractMagnetization

The spontaneous (symmetry-broken) order-parameter magnitude `M(T)` in
the ordered phase; identically zero above `T_c`.  A scalar — the
magnitude, not a spin component — so `tensor_rank == 0`; its critical
exponent is β (see `critical_scaling`).
"""
struct SpontaneousMagnetization <: AbstractMagnetization end
export SpontaneousMagnetization

function _axistuple(I)
    return if (I isa Tuple && !isempty(I) && all(a -> a isa Symbol, I))
        I
    else
        error(
            "index parameter must be a non-empty Tuple of spin-axis Symbols, got $(repr(I))"
        )
    end
end

"""
    Susceptibility{I}() <: AbstractSusceptibility
    Susceptibility(α, β₁, …, βₙ)          # each a Symbol

Susceptibility of **arbitrary response order** — the `n`-th order term
of the order-parameter response to its conjugate field,

`χ⁽ⁿ⁾_{α; β₁…βₙ} = ∂ⁿ⟨M_α⟩ / ∂h_{β₁}…∂h_{βₙ}`,

a rank-`(n+1)` tensor in [`SpinAxis`](@ref) space whose index parameter
`I = (α, β₁, …, βₙ)` carries one response direction `α` and `n` field
directions.  The response order is `n = length(I) − 1`
([`response_order`](@ref)):

- `Susceptibility(:x, :y)` — **linear** `χ_xy = ∂M_x/∂h_y` (order 1),
  the off-diagonal component the fused `SusceptibilityXX/ZZ` names could
  not express;
- `Susceptibility(:x, :y, :z)` — **second-order nonlinear**
  `χ⁽²⁾_{x;yz} = ∂²M_x/∂h_y∂h_z` (order 2);
- `Susceptibility(:x, :x, :x, :x)` — third-order `χ⁽³⁾`, and so on.

The genealogy is recursive: `χ⁽ⁿ⁾ ⟵ χ⁽ⁿ⁻¹⁾ ⟵ … ⟵ M ⟵ F`
(`derivative_edge`), so `derivative_order(χ⁽ⁿ⁾, MagneticField) == n + 1`.
The linear component's defining identity is [`SusceptibilityFDT`](@ref)
`χ_AB = β·Cov(M_A, M_B)`.
"""
struct Susceptibility{I} <: AbstractSusceptibility
    function Susceptibility{I}() where {I}
        return (
            length(_axistuple(I)) >= 2 || error(
                "Susceptibility needs ≥2 indices (1 response + ≥1 field), got $(repr(I))",
            );
            new{I}()
        )
    end
end
Susceptibility(idx::Symbol...) = Susceptibility{idx}()
tensor_rank(::Type{Susceptibility{I}}) where {I} = length(I)
index_spaces(::Type{Susceptibility{I}}) where {I} = ntuple(_ -> SpinAxis(), length(I))
indices(::Type{Susceptibility{I}}) where {I} = I
response_order(::Type{Susceptibility{I}}) where {I} = length(I) - 1
export Susceptibility

# ─── Two-point correlations (tensors in spin space) ─────────────────────

"""
    SpinCorrelation{A,B}() <: AbstractTwoPointCorrelation
    SpinCorrelation(a::Symbol, b::Symbol)

Two-point spin correlation `⟨S^A_i S^B_j⟩` — a rank-2 tensor in
[`SpinAxis`](@ref) space (`SpinCorrelation(:z, :z)` replaces the old
`ZZCorrelation`).  At criticality its decay is governed by the anomalous
dimension η — see the `correlation_decay` correspondence.
"""
struct SpinCorrelation{A,B} <: AbstractTwoPointCorrelation
    SpinCorrelation{A,B}() where {A,B} = (_axis(A); _axis(B); new{A,B}())
end
SpinCorrelation(a::Symbol, b::Symbol) = SpinCorrelation{a,b}()
tensor_rank(::Type{<:SpinCorrelation}) = 2
index_spaces(::Type{<:SpinCorrelation}) = (SpinAxis(), SpinAxis())
indices(::Type{SpinCorrelation{A,B}}) where {A,B} = (A, B)
export SpinCorrelation

"""
    Conductivity{I}() <: AbstractQuantity
    Conductivity(μ, ν₁, …, νₙ)            # each a Symbol

The **DC (static)** electrical conductivity of **arbitrary response
order** — the `n`-th order current response
`j_μ = Σ σ⁽ⁿ⁾_{μ; ν₁…νₙ} E_{ν₁}…E_{νₙ}`, a rank-`(n+1)` tensor in
[`SpatialDirection`](@ref) space with one current direction `μ` and `n`
field directions.  `response_order = length(I) − 1`:

- `Conductivity(:x, :y)` — **linear** `σ_xy` (order 1); its Hall
  component is quantized by [`TKNN`](@ref);
- `Conductivity(:x, :y, :z)` — **second-order** `σ⁽²⁾`, and so on.

This is the zero-frequency response (`frequency_arguments == 0`), the
current-channel analogue of the static [`Susceptibility`](@ref); like it,
it carries intrinsic permutation symmetry over its field indices (at zero
frequency).  Its `ω → 0` limit fixes it from the frequency-resolved AC
[`DynamicalConductivity`](@ref) `σ⁽ⁿ⁾(ω₁, …, ωₙ)` (optical `σ(ω)`, the
photogalvanic `σ⁽²⁾(ω₁, ω₂)`, Drude / f-sum rule) — the current-channel
mirror of [`DynamicalSusceptibility`](@ref).
"""
struct Conductivity{I} <: AbstractQuantity
    function Conductivity{I}() where {I}
        return (
            length(_axistuple(I)) >= 2 || error(
                "Conductivity needs ≥2 indices (1 current + ≥1 field), got $(repr(I))"
            );
            new{I}()
        )
    end
end
Conductivity(idx::Symbol...) = Conductivity{idx}()
tensor_rank(::Type{Conductivity{I}}) where {I} = length(I)
index_spaces(::Type{Conductivity{I}}) where {I} = ntuple(_ -> SpatialDirection(), length(I))
indices(::Type{Conductivity{I}}) where {I} = I
response_order(::Type{Conductivity{I}}) where {I} = length(I) - 1
export Conductivity

"""
    DynamicalConductivity{I}() <: AbstractQuantity
    DynamicalConductivity(μ, ν₁, …, νₙ)           # each a Symbol

The **AC (frequency-resolved)** electrical conductivity of arbitrary
response order — the current-channel mirror of
[`DynamicalSusceptibility`](@ref) and the frequency-resolved counterpart
of the DC [`Conductivity`](@ref) (its `ω → 0` limit).

The linear `DynamicalConductivity(:x, :y)` is the optical conductivity
`σ_xy(ω)` (Drude peak, f-sum rule, Kramers–Kronig between Re and Im).
The `n`-th order term `DynamicalConductivity(μ, ν₁, …, νₙ)` is
`σ⁽ⁿ⁾_{μ; ν₁…νₙ}(ω₁, …, ωₙ)`: the field acts at `n` distinct times, so
the response is intrinsically **multi-time** — `frequency_arguments ==
n == response_order`.  `DynamicalConductivity(:x, :y, :z)` is the
second-order `σ⁽²⁾(ω₁, ω₂)` of the photogalvanic / second-harmonic
response.  Its microscopic Kubo expression is the retarded `n`-time
current–current correlation ([`CurrentCorrelation`](@ref); see
`structure/spectral.jl`).
"""
struct DynamicalConductivity{I} <: AbstractQuantity
    function DynamicalConductivity{I}() where {I}
        return (
            length(_axistuple(I)) >= 2 || error(
                "DynamicalConductivity needs ≥2 indices (1 current + ≥1 field), got $(repr(I))",
            );
            new{I}()
        )
    end
end
DynamicalConductivity(idx::Symbol...) = DynamicalConductivity{idx}()
tensor_rank(::Type{DynamicalConductivity{I}}) where {I} = length(I)
function index_spaces(::Type{DynamicalConductivity{I}}) where {I}
    return ntuple(_ -> SpatialDirection(), length(I))
end
indices(::Type{DynamicalConductivity{I}}) where {I} = I
response_order(::Type{DynamicalConductivity{I}}) where {I} = length(I) - 1
# multi-time: an n-th order AC response depends on n frequencies
frequency_arguments(::Type{DynamicalConductivity{I}}) where {I} = length(I) - 1
export DynamicalConductivity

"""
    CurrentCorrelation{I}() <: AbstractQuantity
    CurrentCorrelation(μ, ν₁, …, νₙ)              # each a Symbol

The `n`-time current–current correlation — the microscopic Kubo kernel of
the [`DynamicalConductivity`](@ref), the current-channel analogue of the
[`DynamicalCorrelation`](@ref).  The linear `CurrentCorrelation(:x, :y)`
is the two-point `⟨j_x(t) j_y(0)⟩` whose retarded part gives `σ_xy(ω)`;
the order-`n` term is the `(n+1)`-point current correlation with `n`
independent time differences (`frequency_arguments == n`), matching the
order of the conductivity it feeds (order-faithful Kubo edge).
"""
struct CurrentCorrelation{I} <: AbstractQuantity
    function CurrentCorrelation{I}() where {I}
        return (
            length(_axistuple(I)) >= 2 || error(
                "CurrentCorrelation needs ≥2 operators (1 + n for order n), got $(repr(I))",
            );
            new{I}()
        )
    end
end
CurrentCorrelation(idx::Symbol...) = CurrentCorrelation{idx}()
tensor_rank(::Type{CurrentCorrelation{I}}) where {I} = length(I)
function index_spaces(::Type{CurrentCorrelation{I}}) where {I}
    return ntuple(_ -> SpatialDirection(), length(I))
end
indices(::Type{CurrentCorrelation{I}}) where {I} = I
response_order(::Type{CurrentCorrelation{I}}) where {I} = length(I) - 1
frequency_arguments(::Type{CurrentCorrelation{I}}) where {I} = length(I) - 1
export CurrentCorrelation

"""
    CurrentNoise{I}() <: AbstractQuantity
    CurrentNoise(μ, ν)                            # each a Symbol

The (symmetrized) current-noise spectral density `S^j_μν(q, ω)` — the
current-channel structure factor: the space-time Fourier transform of the
[`CurrentCorrelation`](@ref) (mirroring
[`DynamicalStructureFactor`](@ref) ↔ [`DynamicalCorrelation`](@ref)) and
the **fluctuation** partner of the dissipative `Re σ_μν(ω)` via the
Johnson–Nyquist fluctuation–dissipation theorem (Nyquist, Phys. Rev. 32,
110 (1928); Callen & Welton, Phys. Rev. 83, 34 (1951)).  `frequency_arguments == 1`.
"""
struct CurrentNoise{I} <: AbstractQuantity
    function CurrentNoise{I}() where {I}
        return (
            length(_axistuple(I)) == 2 ||
                error("CurrentNoise is a rank-2 tensor S^j_μν (2 indices), got $(repr(I))");
            new{I}()
        )
    end
end
CurrentNoise(idx::Symbol...) = CurrentNoise{idx}()
tensor_rank(::Type{CurrentNoise{I}}) where {I} = length(I)
index_spaces(::Type{CurrentNoise{I}}) where {I} = ntuple(_ -> SpatialDirection(), length(I))
indices(::Type{CurrentNoise{I}}) where {I} = I
frequency_arguments(::Type{CurrentNoise{I}}) where {I} = 1
export CurrentNoise

# ─── Transport family ───────────────────────────────────────────────────
#
# The linear-transport quantities: the currents that respond to an
# electric field / temperature gradient, and the rank-2 transport
# coefficients that relate them.  The identities linking them
# (Wiedemann–Franz, Mott, Kelvin, Onsager, the optical sum rule) live in
# `relations/transport.jl`.

"""
    ElectricCurrent() <: AbstractQuantity

The electric (charge) current density `j_μ` — a rank-1 vector in
[`SpatialDirection`](@ref) space; the response half of the
[`Conductivity`](@ref) (`j_μ = σ_μν E_ν`) and one of the two currents of
the Onsager transport matrix (with [`HeatCurrent`](@ref)).
"""
struct ElectricCurrent <: AbstractQuantity end
export ElectricCurrent
tensor_rank(::Type{ElectricCurrent}) = 1
index_spaces(::Type{ElectricCurrent}) = (SpatialDirection(),)

"""
    HeatCurrent() <: AbstractQuantity

The heat (thermal energy) current density `j^Q_μ` — a rank-1 vector in
[`SpatialDirection`](@ref) space; the current driven by a temperature
gradient (`j^Q_μ = −κ_μν ∂_ν T` at zero electric current) and the Onsager
partner of the [`ElectricCurrent`](@ref).
"""
struct HeatCurrent <: AbstractQuantity end
export HeatCurrent
tensor_rank(::Type{HeatCurrent}) = 1
index_spaces(::Type{HeatCurrent}) = (SpatialDirection(),)

"""
    DrudeWeight{I}() <: AbstractQuantity
    DrudeWeight(μ, ν)                             # each a Symbol

The Drude weight (charge stiffness) tensor `D_μν` — the coefficient of
the zero-frequency delta in the real optical conductivity,
`Re σ_μν(ω) = π D_μν δ(ω) + σ^reg_μν(ω)` (Scalapino, White & Zhang, Phys.
Rev. B 47, 7995 (1993)).  A rank-2 tensor in [`SpatialDirection`](@ref)
space; `D_μν > 0` signals a (perfect) conductor.  Fixed by the
[`DynamicalConductivity`](@ref) via the optical sum rule.
"""
struct DrudeWeight{I} <: AbstractQuantity
    function DrudeWeight{I}() where {I}
        return (
            length(_axistuple(I)) == 2 ||
                error("DrudeWeight is a rank-2 tensor D_μν (2 indices), got $(repr(I))");
            new{I}()
        )
    end
end
DrudeWeight(idx::Symbol...) = DrudeWeight{idx}()
tensor_rank(::Type{DrudeWeight{I}}) where {I} = length(I)
index_spaces(::Type{DrudeWeight{I}}) where {I} = ntuple(_ -> SpatialDirection(), length(I))
indices(::Type{DrudeWeight{I}}) where {I} = I
export DrudeWeight

"""
    ThermalConductivity{I}() <: AbstractQuantity
    ThermalConductivity(μ, ν)                     # each a Symbol

The (DC) thermal conductivity tensor `κ_μν` — the heat-current response to
a temperature gradient, `j^Q_μ = −κ_μν ∂_ν T`.  Rank-2 in
[`SpatialDirection`](@ref) space; its ratio to the electrical
[`Conductivity`](@ref) is fixed by the Wiedemann–Franz law.
"""
struct ThermalConductivity{I} <: AbstractQuantity
    function ThermalConductivity{I}() where {I}
        return (
            length(_axistuple(I)) == 2 || error(
                "ThermalConductivity is a rank-2 tensor κ_μν (2 indices), got $(repr(I))",
            );
            new{I}()
        )
    end
end
ThermalConductivity(idx::Symbol...) = ThermalConductivity{idx}()
tensor_rank(::Type{ThermalConductivity{I}}) where {I} = length(I)
function index_spaces(::Type{ThermalConductivity{I}}) where {I}
    return ntuple(_ -> SpatialDirection(), length(I))
end
indices(::Type{ThermalConductivity{I}}) where {I} = I
export ThermalConductivity

"""
    Thermopower{I}() <: AbstractQuantity
    Thermopower(μ, ν)                             # each a Symbol

The thermopower (Seebeck coefficient) tensor `S_μν` — the electric field
generated per unit temperature gradient at zero current,
`E_μ = S_μν ∂_ν T`.  Rank-2 in [`SpatialDirection`](@ref) space; fixed by
the Mott formula and linked to the [`PeltierCoefficient`](@ref) by the
Kelvin relation.
"""
struct Thermopower{I} <: AbstractQuantity
    function Thermopower{I}() where {I}
        return (
            length(_axistuple(I)) == 2 ||
                error("Thermopower is a rank-2 tensor S_μν (2 indices), got $(repr(I))");
            new{I}()
        )
    end
end
Thermopower(idx::Symbol...) = Thermopower{idx}()
tensor_rank(::Type{Thermopower{I}}) where {I} = length(I)
index_spaces(::Type{Thermopower{I}}) where {I} = ntuple(_ -> SpatialDirection(), length(I))
indices(::Type{Thermopower{I}}) where {I} = I
export Thermopower

"""
    PeltierCoefficient{I}() <: AbstractQuantity
    PeltierCoefficient(μ, ν)                      # each a Symbol

The Peltier coefficient tensor `Π_μν` — the heat current carried per unit
electric current, `j^Q_μ = Π_μν j_ν`.  Rank-2 in [`SpatialDirection`](@ref)
space; the Kelvin (second Thomson) relation ties it to the
[`Thermopower`](@ref), `Π = T S`.
"""
struct PeltierCoefficient{I} <: AbstractQuantity
    function PeltierCoefficient{I}() where {I}
        return (
            length(_axistuple(I)) == 2 || error(
                "PeltierCoefficient is a rank-2 tensor Π_μν (2 indices), got $(repr(I))"
            );
            new{I}()
        )
    end
end
PeltierCoefficient(idx::Symbol...) = PeltierCoefficient{idx}()
tensor_rank(::Type{PeltierCoefficient{I}}) where {I} = length(I)
function index_spaces(::Type{PeltierCoefficient{I}}) where {I}
    return ntuple(_ -> SpatialDirection(), length(I))
end
indices(::Type{PeltierCoefficient{I}}) where {I} = I
export PeltierCoefficient

# ─── Dynamical & spectral quantities ────────────────────────────────────
#
# The frequency-resolved family.  These tags name the quantities; the
# identities relating them (`A = −Im G^R/π`, Dyson, DOS = BZ-average of A,
# detailed balance, the NMR relaxation relations) live in
# `relations/spectral.jl` and the transform/sum genealogy in
# `structure/spectral.jl`.  Actually *evaluating* the ω-dependence
# (Kramers–Kronig, analytic continuation) is out of scope here — see the
# functional sibling (issue #14).

"""
    AbstractPropagator <: AbstractQuantity

Single-particle propagators — retarded/advanced/Matsubara Green's
functions and the self-energy — the `(q, ω)`-resolved objects the Dyson
equation relates.
"""
abstract type AbstractPropagator <: AbstractQuantity end
export AbstractPropagator

"""
    RetardedGreensFunction() <: AbstractPropagator

The retarded single-particle Green's function `G^R(q, ω)`.  Its spectral
representation `A = −Im G^R/π` and the Dyson equation
`G^{-1} = G₀^{-1} − Σ` are in `relations/spectral.jl`.
"""
struct RetardedGreensFunction <: AbstractPropagator end
export RetardedGreensFunction
# G_ab(q,ω): rank-2 in orbital space
tensor_rank(::Type{RetardedGreensFunction}) = 2
index_spaces(::Type{RetardedGreensFunction}) = (OrbitalIndex(), OrbitalIndex())
frequency_arguments(::Type{RetardedGreensFunction}) = 1

"""
    SelfEnergy() <: AbstractPropagator

The single-particle self-energy `Σ(q, ω)` — the Dyson correction
`G^{-1} = G₀^{-1} − Σ` between the bare and full propagators.
"""
struct SelfEnergy <: AbstractPropagator end
export SelfEnergy
tensor_rank(::Type{SelfEnergy}) = 2
index_spaces(::Type{SelfEnergy}) = (OrbitalIndex(), OrbitalIndex())
frequency_arguments(::Type{SelfEnergy}) = 1

"""
    SpectralFunction() <: AbstractQuantity

The single-particle spectral function `A(q, ω) = −(1/π) Im G^R(q, ω)`,
normalized by `∫ A(q, ω) dω = 1`.
"""
struct SpectralFunction <: AbstractQuantity end
export SpectralFunction
tensor_rank(::Type{SpectralFunction}) = 2
index_spaces(::Type{SpectralFunction}) = (OrbitalIndex(), OrbitalIndex())
frequency_arguments(::Type{SpectralFunction}) = 1

"""
    DensityOfStates() <: AbstractQuantity

The density of states `ρ(ω) = (1/N) Σ_q A(q, ω)` — the Brillouin-zone
average of the [`SpectralFunction`](@ref).
"""
struct DensityOfStates <: AbstractQuantity end
export DensityOfStates
frequency_arguments(::Type{DensityOfStates}) = 1

"""
    DynamicalCorrelation{I}() <: AbstractQuantity
    DynamicalCorrelation(α, β₁, …, βₙ)             # each a Symbol

The space-and-time-resolved correlation of **arbitrary order** — the
microscopic kernel of the Kubo response, carrying the **same order** as
the [`DynamicalSusceptibility`](@ref) it feeds.

The linear `DynamicalCorrelation(:x, :y)` is the two-point
`⟨A^x(r, t) A^y(0, 0)⟩` whose space-time Fourier transform is the
[`DynamicalStructureFactor`](@ref) `S(q, ω)` — one time difference,
`frequency_arguments == 1`.

The `n`-th order term `DynamicalCorrelation(α, β₁, …, βₙ)` is the
`(n+1)`-point function `⟨A^α(t) A^{β₁}(t₁) ⋯ A^{βₙ}(tₙ)⟩` — `n+1`
operators at `n` independent time differences, so it is intrinsically
**n-time** (`frequency_arguments == n == response_order`).  Its `n`-fold
nested-commutator (retarded) part is exactly the Kubo kernel of the
order-`n` `DynamicalSusceptibility(α, β₁, …, βₙ)` (Kubo, J. Phys. Soc.
Jpn. 12, 570 (1957)): an `n`-th order response is an `n`-time
correlation.
"""
struct DynamicalCorrelation{I} <: AbstractQuantity
    function DynamicalCorrelation{I}() where {I}
        return (
            length(_axistuple(I)) >= 2 || error(
                "DynamicalCorrelation needs ≥2 operators (1 + n fields for order n), got $(repr(I))",
            );
            new{I}()
        )
    end
end
DynamicalCorrelation(idx::Symbol...) = DynamicalCorrelation{idx}()
export DynamicalCorrelation
tensor_rank(::Type{DynamicalCorrelation{I}}) where {I} = length(I)
function index_spaces(::Type{DynamicalCorrelation{I}}) where {I}
    return ntuple(_ -> SpinAxis(), length(I))
end
indices(::Type{DynamicalCorrelation{I}}) where {I} = I
response_order(::Type{DynamicalCorrelation{I}}) where {I} = length(I) - 1
# n-time: an (n+1)-point correlation has n independent time differences,
# matching the n frequencies of the order-n dynamical susceptibility.
frequency_arguments(::Type{DynamicalCorrelation{I}}) where {I} = length(I) - 1

"""
    DynamicalStructureFactor() <: AbstractStructureFactor

The dynamical structure factor `S(q, ω)` — the space-time Fourier
transform of the [`DynamicalCorrelation`](@ref); obeys detailed balance
`S(q, −ω) = e^{−βω} S(q, ω)` and the fluctuation–dissipation link to the
[`DynamicalSusceptibility`](@ref).
"""
struct DynamicalStructureFactor <: AbstractStructureFactor end
export DynamicalStructureFactor
# S_αβ(q,ω): rank-2 in spin space
tensor_rank(::Type{DynamicalStructureFactor}) = 2
index_spaces(::Type{DynamicalStructureFactor}) = (SpinAxis(), SpinAxis())
frequency_arguments(::Type{DynamicalStructureFactor}) = 1

"""
    StaticStructureFactor() <: AbstractStructureFactor

The static (equal-time) structure factor `S(q)` — the frequency integral
of the [`DynamicalStructureFactor`](@ref), `S(q) = ∫ S(q, ω) dω/(2π)`
(Van Hove, Phys. Rev. 95, 249 (1954)).  Its `q → 0` limit fixes the
static susceptibility (`χ = β S(q→0)`, classical).  Rank-2 in spin space,
one frequency integrated out (`frequency_arguments == 0`).
"""
struct StaticStructureFactor <: AbstractStructureFactor end
export StaticStructureFactor
tensor_rank(::Type{StaticStructureFactor}) = 2
index_spaces(::Type{StaticStructureFactor}) = (SpinAxis(), SpinAxis())

"""
    DynamicalSusceptibility{I}() <: AbstractSusceptibility
    DynamicalSusceptibility(α, β₁, …, βₙ)          # each a Symbol

The dynamical susceptibility of **arbitrary response order** — the
frequency-domain (multi-time) counterpart of the static
[`Susceptibility`](@ref).  The linear `DynamicalSusceptibility(:x, :y)`
is `χ_xy(ω)`, one frequency argument, and its imaginary part `χ''(q, ω)`
is the dissipative response of the fluctuation–dissipation theorem and
the NMR relaxation rate.

The `n`-th order term `DynamicalSusceptibility(α, β₁, …, βₙ)` is
`χ⁽ⁿ⁾_{α;β₁…βₙ}(ω₁, …, ωₙ)`: the field is applied at `n` distinct times,
so the response is intrinsically **multi-time** — `frequency_arguments
== n` (`response_order`).  `DynamicalSusceptibility(:x, :y, :z)` is the
second-order `χ⁽²⁾(ω₁, ω₂)` of two-dimensional coherent spectroscopy
(Wan & Armitage, Phys. Rev. Lett. 122, 257401 (2019)).  Its microscopic
Kubo expression is the `n`-fold nested-commutator response function
(Kubo, J. Phys. Soc. Jpn. 12, 570 (1957)); see `structure/spectral.jl`.

The static `Susceptibility{I}` of the same order is the zero-frequency
limit, `χ⁽ⁿ⁾(0, …, 0)`.
"""
struct DynamicalSusceptibility{I} <: AbstractSusceptibility
    function DynamicalSusceptibility{I}() where {I}
        return (
            length(_axistuple(I)) >= 2 || error(
                "DynamicalSusceptibility needs ≥2 indices (1 response + ≥1 field), got $(repr(I))",
            );
            new{I}()
        )
    end
end
DynamicalSusceptibility(idx::Symbol...) = DynamicalSusceptibility{idx}()
tensor_rank(::Type{DynamicalSusceptibility{I}}) where {I} = length(I)
function index_spaces(::Type{DynamicalSusceptibility{I}}) where {I}
    return ntuple(_ -> SpinAxis(), length(I))
end
indices(::Type{DynamicalSusceptibility{I}}) where {I} = I
response_order(::Type{DynamicalSusceptibility{I}}) where {I} = length(I) - 1
# multi-time: an n-th order dynamical response depends on n frequencies
frequency_arguments(::Type{DynamicalSusceptibility{I}}) where {I} = length(I) - 1
export DynamicalSusceptibility

"""
    NMRSpinRelaxationRate() <: AbstractQuantity

The NMR spin–lattice relaxation rate `1/T₁` — set by the low-frequency
limit of the dissipative dynamical susceptibility (Moriya),
`1/T₁ ∝ T · lim_{ω→0} Σ_q |A_hf(q)|² χ''(q, ω)/ω`.
"""
struct NMRSpinRelaxationRate <: AbstractQuantity end
export NMRSpinRelaxationRate

"""
    NMRRelaxationExponent() <: AbstractQuantity

The low-temperature scaling exponent `θ_NMR` of `1/T₁ ∝ T^{θ_NMR}`,
fixed by the operator scaling dimension via `θ_NMR = 2Δ_op − 1`.
"""
struct NMRRelaxationExponent <: AbstractQuantity end
export NMRRelaxationExponent

"""
    ScalingDimension() <: AbstractQuantity

The scaling dimension `Δ_op` of a local operator at a quantum critical
point — the input to dynamical scaling relations such as the NMR
exponent `θ_NMR = 2Δ_op − 1`.
"""
struct ScalingDimension <: AbstractQuantity end
export ScalingDimension

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

Correlation length ξ (units of lattice spacing).  In a gapped phase it
is set by the gap and velocity, `ξ = v/Δ` ([`CorrelationLengthGap`](@ref)).
"""
struct CorrelationLength <: AbstractQuantity end
export CorrelationLength

"""
    MassGap() <: AbstractGap

The spectral (mass) gap `Δ = E₁ − E₀` between the ground state and the
first excitation.  Sets the correlation length `ξ = v/Δ` in a gapped
phase, and vanishes as `Δ ∼ ξ^{−z}` (dynamical exponent `z`) on approach
to a quantum critical point.
"""
struct MassGap <: AbstractGap end
export MassGap

"""
    DynamicalExponent() <: AbstractQuantity

The dynamical critical exponent `z` relating spatial and temporal
scaling at a quantum critical point, `Δ ∼ ξ^{−z}` (equivalently
`ω ∼ k^z`).  `z = 1` for a Lorentz-invariant (relativistic) critical
point.
"""
struct DynamicalExponent <: AbstractQuantity end
export DynamicalExponent

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

# ─── Entanglement ────────────────────────────────────────────────────────

"""
    VonNeumannEntropy() <: AbstractEntanglementMeasure

The von Neumann entanglement entropy `S = −Tr(ρ_A ln ρ_A)` of a
subsystem — the `n → 1` limit of the [`RenyiEntropy`](@ref).  In a gapped
phase it obeys an area law (Eisert, Cramer & Plenio, Rev. Mod. Phys. 82,
277 (2010)); at a 1D critical point it grows logarithmically with the
subsystem size, `S = (c/3) ln ℓ` (Calabrese & Cardy, J. Stat. Mech.
(2004) P06002).
"""
struct VonNeumannEntropy <: AbstractEntanglementMeasure end
export VonNeumannEntropy

"""
    RenyiEntropy() <: AbstractEntanglementMeasure

The Rényi entanglement entropy `S_n = (1−n)⁻¹ ln Tr(ρ_A^n)`.  The order
`n` is supplied at use; `n = 2` is fixed by the [`Purity`](@ref)
(`S_2 = −ln Tr ρ_A²`), and `n → 1` recovers the [`VonNeumannEntropy`](@ref).
"""
struct RenyiEntropy <: AbstractEntanglementMeasure end
export RenyiEntropy

"""
    Purity() <: AbstractQuantity

The purity `Tr(ρ_A²) ∈ (0, 1]` of a (reduced) density matrix — `1` for a
pure state, `1/d` for the maximally mixed one.  Fixes the Rényi-2 entropy
via `S_2 = −ln Tr ρ_A²`.
"""
struct Purity <: AbstractQuantity end
export Purity

"""
    ChernNumber() <: AbstractQuantity

The (first) Chern number `C ∈ ℤ` of a set of bands — the Brillouin-zone
integral of the Berry curvature, `C = (1/2π) ∫_BZ Ω(k) d²k` (Thouless,
Kohmoto, Nightingale & den Nijs, Phys. Rev. Lett. 49, 405 (1982)).  It
sets the quantized Hall conductance ([`TKNN`](@ref)) and, via the
bulk–boundary correspondence, the number of chiral edge modes.
"""
struct ChernNumber <: AbstractQuantity end
export ChernNumber

"""
    BerryCurvature() <: AbstractQuantity

The Berry curvature `Ω(k)` of a band — the momentum-space field strength
`Ω = ∂_{k_x} A_y − ∂_{k_y} A_x` of the Berry connection (Berry, Proc. R.
Soc. Lond. A 392, 45 (1984)).  Its Brillouin-zone integral is the
[`ChernNumber`](@ref); it also drives the intrinsic anomalous Hall
effect (Xiao, Chang & Niu, Rev. Mod. Phys. 82, 1959 (2010)).

Note (scope): the Berry curvature is the *imaginary* part of the quantum
geometric tensor; the real part (the quantum metric) and the mixed-state
/ Uhlmann generalizations are deliberately out of this package's scope.
"""
struct BerryCurvature <: AbstractQuantity end
export BerryCurvature

"""
    BoundaryModeCount() <: AbstractQuantity

The number of protected boundary (edge / surface) modes of a
topological phase — fixed by the bulk topological invariant through the
bulk–boundary correspondence, `n = |ν|` (Hasan & Kane, Rev. Mod. Phys.
82, 3045 (2010)).  See [`BulkBoundary`](@ref).
"""
struct BoundaryModeCount <: AbstractQuantity end
export BoundaryModeCount
