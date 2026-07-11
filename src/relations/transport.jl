# relations/transport.jl — the linear-transport identities.
#
# Model-independent relations among the transport coefficients
# (electrical / thermal conductivity, thermopower, Peltier, Drude weight):
# Wiedemann–Franz, the Mott thermopower formula, the Kelvin (second
# Thomson) relation, Onsager reciprocity, and the optical (Drude/regular)
# sum rule.  Domain tag :transport throughout.
#
# Natural units (k_B = e = 1) as elsewhere in the package; the Sommerfeld
# constants are carried as supplied values (e.g. the Lorenz number `L0`)
# so the caller checks a computed coefficient against the identity.
#
# References (doiget-verified, docs/references.bib): Onsager, Phys. Rev.
# 37, 405 (1931) and 38, 2265 (1931); Cutler & Mott, Phys. Rev. 181, 1336
# (1969); Scalapino, White & Zhang, Phys. Rev. B 47, 7995 (1993).

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
