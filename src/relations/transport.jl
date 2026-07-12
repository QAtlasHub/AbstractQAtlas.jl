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
# (1969); Scalapino, White & Zhang, Phys. Rev. B 47, 7995 (1993); Drude,
# Ann. Phys. 306, 566 (1900); Einstein, Ann. Phys. 322, 549 (1905); Hall,
# Am. J. Math. 2, 287 (1879).

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
@relation :transport WiedemannFranz(κ, σ, T, L0) = κ - L0 * σ * T

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
@relation :transport MottFormula(S, dlnσ_dε, T) = S + (π^2 / 3) * T * dlnσ_dε

"""
    KelvinRelation <: AbstractRelation

The Kelvin (second Thomson) relation — a consequence of Onsager
reciprocity (Onsager, Phys. Rev. 37, 405 (1931)) — tying the Peltier
coefficient to the thermopower,

`Π = T · S`.

Variables: `Π`, `S`, `T`.
"""
@relation :transport KelvinRelation(Π, S, T) = Π - T * S

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
@relation :transport OpticalSumRule(sigma_integral, D, W_reg) =
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
@relation :transport CurrentNoiseFDT(S_j, Reσ, ω, β) = S_j - ω * coth(β * ω / 2) * Reσ

"""
    MobilityConductivity <: AbstractRelation

The Drude relation between conductivity, carrier density and mobility,

`σ = n e μ`,

the current-channel definition of the [`Mobility`](@ref) `μ`
([`CarrierDensity`](@ref) `n`, charge `e`).

Variables: `σ`, `n`, `e`, `μ`.
"""
@relation :transport MobilityConductivity(σ, n, e, μ) = σ - n * e * μ

"""
    DrudeMobility <: AbstractRelation

The Drude mobility from the transport scattering time and effective mass
(Drude, Ann. Phys. 306, 566 (1900)),

`μ = e τ / m`,

([`ScatteringTime`](@ref) `τ`, [`EffectiveMass`](@ref) `m`).  Composed
with [`MobilityConductivity`](@ref) it gives `σ = n e² τ / m`.

Variables: `μ`, `e`, `τ`, `m`.
"""
@relation :transport DrudeMobility(μ, e, τ, m) = μ - e * τ / m

"""
    EinsteinRelation <: AbstractRelation

The Einstein (Einstein–Smoluchowski) relation between mobility and the
diffusion constant (Einstein, Ann. Phys. 322, 549 (1905)),

`μ = e D / k_B T = e D β`,

the fluctuation–dissipation link for transport ([`DiffusionConstant`](@ref)
`D`).  With [`DrudeMobility`](@ref) it fixes `D = k_B T τ / m`.

Variables: `μ`, `e`, `D`, and `β` (or `T`).
"""
@relation :transport EinsteinRelation(μ, e, D, β) = μ - e * D * β

"""
    SingleBandHall <: AbstractRelation

The single-band Hall coefficient (Hall, Am. J. Math. 2, 287 (1879)),

`R_H = 1 / (n e)`,

fixing the [`HallCoefficient`](@ref) `R_H` (and its sign) from the carrier
density and charge (`R_H·n·e = 1`; the sign of `e` gives the carrier sign).

Variables: `R_H`, `n`, `e`.
"""
@relation :transport SingleBandHall(R_H, n, e) = R_H * n * e - 1

"""
    HallAngle <: AbstractRelation

The Hall angle from the conductivity tensor,

`tan θ_H = σ_xy / σ_xx`  (`= ω_c τ` in the Drude picture),

the ratio of the transverse (Hall) to longitudinal conductivity.

Variables: `tanθ_H`, `σxy`, `σxx`.
"""
@relation :transport HallAngle(tanθ_H, σxy, σxx) = tanθ_H - σxy / σxx
