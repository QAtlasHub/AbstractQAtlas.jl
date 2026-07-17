# relations/transport.jl — the linear-transport identities.
#
# Model-independent relations among the transport coefficients
# (electrical / thermal conductivity, thermopower, Peltier, Drude weight):
# Wiedemann–Franz, the Mott thermopower formula, the Kelvin (second
# Thomson) relation, Onsager reciprocity, and the optical (Drude/regular)
# sum rule.  Domain tag :transport throughout.
#
# Natural units (k_B = 1) as elsewhere in the package; Sommerfeld
# constants are carried as supplied values (e.g. the Lorenz number `L0`).
# The elementary charge `e` is kept as an explicit variable in the
# carrier / Hall / Drude relations, where it clarifies the identity (pass
# `e = 1` for natural units).
#
# References (doiget-verified, docs/references.bib): Onsager, Phys. Rev.
# 37, 405 (1931) and 38, 2265 (1931); Cutler & Mott, Phys. Rev. 181, 1336
# (1969); Scalapino, White & Zhang, Phys. Rev. B 47, 7995 (1993); Einstein,
# Ann. Phys. 322, 549 (1905).
#
# NOTE: model-SPECIFIC transport relations (the Drude mobility μ=eτ/m, the
# single-band Hall coefficient R_H=1/ne) are deliberately NOT here — they
# assume a particular model / band structure, so they belong to the
# implementing atlas (QAtlas), not this universal layer.

"""
    WiedemannFranz <: AbstractRelation

The Wiedemann–Franz law: the ratio of the thermal to the electrical
conductivity is the temperature times the Lorenz number,

`κ = L₀ · σ · T`,

with the Sommerfeld value `L₀ = π²/3` (in units `k_B = e = 1`; i.e.
`π²k_B²/3e²`).  A diagonal-component statement (`κ_xx`, `σ_xx`); the ratio
`κ/(σT)` is the caller-supplied Lorenz number `L0`, checked against the
Sommerfeld constant.

Variables: `κ`, `σ`, `T`, `L0`.
"""
@relation :transport WiedemannFranz(
    κ::ThermalConductivity{(:x, :x)}, σ::Conductivity{(:x, :x)}, T::Temperature, L0
) = κ - L0 * σ * T

"""
    MottFormula <: AbstractRelation

The Mott formula for the diffusive thermopower (Cutler & Mott, Phys. Rev.
181, 1336 (1969)),

`S = −(π²/3) · T · d ln σ(ε)/dε |_{ε_F}`,

(natural units `k_B = e = 1`, electron-like carriers), fixing the Seebeck
coefficient from the energy derivative of the conductivity at the Fermi
level.

Variables: `S`, `dlnσ_dε` = `d ln σ/dε |_{ε_F}`, `T`.
"""
@relation :transport MottFormula(S::Thermopower{(:x, :x)}, dlnσ_dε, T::Temperature) =
    S + (π^2 / 3) * T * dlnσ_dε

"""
    KelvinRelation <: AbstractRelation

The Kelvin (second Thomson) relation — a consequence of Onsager
reciprocity (Onsager, Phys. Rev. 37, 405 (1931)) — tying the Peltier
coefficient to the thermopower,

`Π = T · S`.

Variables: `Π`, `S`, `T`.
"""
@relation :transport KelvinRelation(
    Π::PeltierCoefficient{(:x, :x)}, S::Thermopower{(:x, :x)}, T::Temperature
) = Π - T * S

"""
    OnsagerReciprocity <: AbstractRelation

Onsager reciprocity (Onsager, Phys. Rev. 37, 405 (1931); 38, 2265
(1931)): in the absence of a magnetic field the linear-transport matrix is
symmetric,

`L_{μν} = L_{νμ}`,

so a transport coefficient equals its index-transposed partner (`σ_xy =
σ_yx`, `κ_xy = κ_yx`, …).  In a field `B` the Onsager–Casimir form
`L_{μν}(B) = L_{νμ}(−B)` holds instead.

Variables: `L_μν`, `L_νμ`.
"""
@relation :transport OnsagerReciprocity(L_μν, L_νμ) = L_μν - L_νμ

"""
    OpticalSumRule <: AbstractRelation

The optical (f-sum) rule decomposing the frequency-integrated real
conductivity into its Drude and regular parts,

`∫ Re σ(ω) dω = π D + W_reg`,

with `D` the [`DrudeWeight`](@ref) (the `δ(ω)` coefficient,
`Re σ = π D δ(ω) + σ^reg`) and `W_reg = ∫ σ^reg(ω) dω` the regular
spectral weight (Scalapino, White & Zhang, Phys. Rev. B 47, 7995 (1993)).
The total `sigma_integral` is the caller-supplied f-sum weight (e.g.
`π n e²/m`, or `−π e²⟨T_kin⟩` on a lattice).

Variables: `sigma_integral`, `D`, `W_reg`.
"""
@relation :transport OpticalSumRule(sigma_integral, D::DrudeWeight{(:x, :x)}, W_reg) =
    sigma_integral - (π * D + W_reg)

"""
    CurrentNoiseFDT <: AbstractRelation

The Johnson–Nyquist fluctuation–dissipation theorem — the **fluctuation**
partner of the dissipative conductivity, the current-channel analogue of
[`DynamicalFDT`](@ref) (`S ↔ χ''`).  The symmetrized current-noise
spectral density ([`CurrentNoise`](@ref)) is fixed by the real
conductivity,

`S^j(ω) = ω · coth(βω/2) · Re σ(ω)`,

(natural units `ℏ = k_B = 1`; Nyquist, Phys. Rev. 32, 110 (1928); Callen &
Welton, Phys. Rev. 83, 34 (1951)).  The classical limit `βω ≪ 1` gives
the white Nyquist noise `S^j = 2 T Re σ` (`ω coth(βω/2) → 2/β`).

Variables: `S_j` = `S^j(ω)`, `Reσ` = `Re σ(ω)`, `ω`, `β` (or `T`).
"""
@relation :transport CurrentNoiseFDT(
    S_j::CurrentNoise{(:x, :x)}, Reσ, ω, β::InverseTemperature
) = S_j - ω * coth(β * ω / 2) * Reσ

"""
    MobilityConductivity <: AbstractRelation

The Drude relation between conductivity, carrier density and mobility,

`σ = n e μ`,

the current-channel definition of the [`Mobility`](@ref) `μ`
([`CarrierDensity`](@ref) `n`, charge `e`).

Variables: `σ`, `n`, `e`, `μ`.
"""
@relation :transport MobilityConductivity(
    σ::Conductivity{(:x, :x)}, n::CarrierDensity, e, μ::Mobility
) = σ - n * e * μ

"""
    EinsteinRelation <: AbstractRelation

The Einstein (Einstein–Smoluchowski) relation between mobility and the
diffusion constant (Einstein, Ann. Phys. 322, 549 (1905)),

`μ = e D / k_B T = e D β`,

the universal fluctuation–dissipation link for transport
([`DiffusionConstant`](@ref) `D`).

Variables: `μ`, `e`, `D`, and `β` (or `T`).
"""
@relation :transport EinsteinRelation(
    μ::Mobility, e, D::DiffusionConstant, β::InverseTemperature
) = μ - e * D * β

"""
    HallAngle <: AbstractRelation

The Hall angle from the conductivity tensor,

`tan θ_H = σ_xy / σ_xx`  (`= ω_c τ` in the Drude picture),

the ratio of the transverse (Hall) to longitudinal conductivity.

Variables: `tanθ_H`, `σxy`, `σxx`.
"""
@relation :transport HallAngle(
    tanθ_H, σxy::Conductivity{(:x, :y)}, σxx::Conductivity{(:x, :x)}
) = tanθ_H - σxy / σxx

"""
    LongitudinalResistivity <: AbstractRelation

The longitudinal resistivity from the 2×2 magnetotransport tensor
inversion `ρ = σ⁻¹`,

`ρ_xx = σ_xx / (σ_xx² + σ_xy²)`,

(convention-free — the diagonal element of the inverse).  In a
dissipationless Hall state (`σ_xx = 0`) it vanishes.

Variables: `ρxx`, `σxx`, `σxy`.
"""
@relation :transport LongitudinalResistivity(
    ρxx::Resistivity{(:x, :x)}, σxx::Conductivity{(:x, :x)}, σxy::Conductivity{(:x, :y)}
) = ρxx - σxx / (σxx^2 + σxy^2)

"""
    HallResistivity <: AbstractRelation

The Hall resistivity from the 2×2 magnetotransport tensor inversion,

`ρ_xy = σ_xy / (σ_xx² + σ_xy²)`,

(standard quantum-Hall sign convention, `ρ_xy` and `σ_xy` carrying the
same sign).  A dissipationless Hall state (`σ_xx = 0`) gives the inverse
Hall conductivity `ρ_xy = 1/σ_xy`.

Variables: `ρxy`, `σxx`, `σxy`.
"""
@relation :transport HallResistivity(
    ρxy::Resistivity{(:x, :y)}, σxx::Conductivity{(:x, :x)}, σxy::Conductivity{(:x, :y)}
) = ρxy - σxy / (σxx^2 + σxy^2)

"""
    CyclotronFrequency <: AbstractRelation

The cyclotron frequency of a carrier in a magnetic field,

`ω_c = e B / m`,

([`MagneticFluxDensity`](@ref) `B`, [`EffectiveMass`](@ref) `m`); sets the
Hall angle `tan θ_H = ω_c τ` ([`HallAngle`](@ref)).

Variables: `ωc`, `e`, `B`, `m`.
"""
@relation :transport CyclotronFrequency(ωc, e, B::MagneticFluxDensity, m::EffectiveMass) =
    ωc - e * B / m

"""
    RighiLeduc <: AbstractRelation

The Righi–Leduc (thermal Hall) effect: the thermal and electrical Hall
conductivities obey the Wiedemann–Franz law in the transverse channel,

`κ_xy = L₀ · T · σ_xy`,

the off-diagonal companion of [`WiedemannFranz`](@ref) (`L₀ = π²/3`).

Variables: `κxy`, `L0`, `T`, `σxy`.
"""
@relation :transport RighiLeduc(
    κxy::ThermalConductivity{(:x, :y)}, L0, T::Temperature, σxy::Conductivity{(:x, :y)}
) = κxy - L0 * T * σxy

"""
    VonKlitzing <: AbstractRelation

The quantized Hall resistance of the integer quantum Hall effect (von
Klitzing, Dorda & Pepper, Phys. Rev. Lett. 45, 494 (1980)),

`R_xy = h / (ν e²) = R_K / ν`,

with the von Klitzing constant `R_K = h/e²` and the integer
[`FillingFactor`](@ref) `ν` (`R_xy·ν·e² = h`).

Variables: `Rxy`, `ν`, `e`, `h`.
"""
@relation :transport VonKlitzing(Rxy::Resistivity{(:x, :y)}, ν::FillingFactor, e, h) =
    Rxy * ν * e^2 - h

"""
    ThermoelectricFigureOfMerit <: AbstractRelation

The dimensionless thermoelectric figure of merit,

`ZT = S² σ T / κ`,

setting the Carnot-fraction efficiency of a thermoelectric
([`Thermopower`](@ref) `S`, electrical [`Conductivity`](@ref) `σ`, total
[`ThermalConductivity`](@ref) `κ`).

Variables: `ZT`, `S`, `σ`, `T`, `κ`.
"""
@relation :transport ThermoelectricFigureOfMerit(
    ZT,
    S::Thermopower{(:x, :x)},
    σ::Conductivity{(:x, :x)},
    T::Temperature,
    κ::ThermalConductivity{(:x, :x)},
) = ZT - S^2 * σ * T / κ

"""
    PowerFactor <: AbstractRelation

The thermoelectric power factor,

`PF = S² σ`,

the numerator of `ZT`; the material's electrical thermoelectric quality
([`Thermopower`](@ref) `S`, [`Conductivity`](@ref) `σ`).

Variables: `PF`, `S`, `σ`.
"""
@relation :transport PowerFactor(PF, S::Thermopower{(:x, :x)}, σ::Conductivity{(:x, :x)}) =
    PF - S^2 * σ

"""
    NernstCoefficient <: AbstractRelation

The Nernst coefficient — the transverse (off-diagonal) thermopower per
unit magnetic field,

`ν_N = S_xy / B`,

the thermoelectric analogue of the Hall effect ([`MagneticFluxDensity`](@ref)
`B`).

Variables: `νN`, `S_xy`, `B`.
"""
@relation :transport NernstCoefficient(
    νN, S_xy::Thermopower{(:x, :y)}, B::MagneticFluxDensity
) = νN - S_xy / B

"""
    IoffeRegel <: AbstractInequality

The Mott–Ioffe–Regel criterion for coherent (metallic) transport
(Ioffe & Regel, Prog. Semicond. 4, 237 (1960)): the mean free path must
exceed the inverse Fermi wavevector,

`k_F ℓ ≥ 1`   (slack `k_F ℓ − 1`).

Its saturation `k_F ℓ ≈ 1` marks the Mott–Ioffe–Regel limit, the breakdown
of Boltzmann quasiparticle transport (the "bad-metal" regime).

Variables: `kFℓ` = `k_F ℓ`.
"""
@inequality :transport IoffeRegel(kFℓ) = kFℓ - 1
