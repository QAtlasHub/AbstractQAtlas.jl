# relations/keldysh.jl — the Keldysh (real-time contour) Green's-function
# structure and its equilibrium fluctuation–dissipation constraint.
#
# On the closed-time contour every two-point function is a 2×2 matrix; the
# retarded–advanced–Keldysh (RAK / Larkin–Ovchinnikov) rotation makes its
# content three independent components `(G^R, G^A, G^K)` built from the
# greater/lesser correlators `G^≷`.  Two of the relations here are ALGEBRAIC
# IDENTITIES that hold in AND out of equilibrium (definitions of the RAK
# components); the third is the CAUSALITY/adjoint tie `G^A = (G^R)†`; and the
# fourth is the physics: in thermal equilibrium the Keldysh component is not
# independent but locked to the spectral part by the FLUCTUATION–DISSIPATION
# theorem,
#
#     G^K(ω) = h(ω) · (G^R(ω) − G^A(ω)),   h(ω) = coth(βω/2) | tanh(βω/2),
#
# the distribution function `h` (bosonic | fermionic) supplied by
# `keldysh_distribution` (structure/keldysh.jl).  The `G^≷` KMS relation
# `G^<(ω) = ζ e^{−βω} G^>(ω)` is the detailed-balance root from which that `h`
# follows (see the test): FDT is a CONSEQUENCE of KMS + the RAK identities,
# not an extra axiom.  This is the real-time sibling of the Matsubara Dyson
# layer in `spectral.jl`; the spectral function `A = −Im G^R/π` bridges them
# ([`SpectralFromKeldysh`](@ref)).
#
# Reference: Keldysh, Sov. Phys. JETP 20, 1018 (1965); the FDT is
# Callen–Welton, [CallenWelton1951](@cite) (already cited).

"""
    KeldyshComponent <: AbstractRelation

Definition of the Keldysh component from the greater/lesser correlators,
`G^K = G^> + G^<` — an identity on the whole contour (equilibrium or not).

Variables: `GK`, `Ggtr` (`G^>`), `Gles` (`G^<`).
"""
@relation :keldysh KeldyshComponent(
    GK::KeldyshGreensFunction, Ggtr::GreaterGreensFunction, Gles::LesserGreensFunction
) = GK - (Ggtr + Gles)

"""
    KeldyshCausality <: AbstractRelation

The retarded–advanced difference equals the greater–lesser difference,
`G^R − G^A = G^> − G^<` — the (un-normalized) spectral weight, an identity
independent of the state.

Variables: `GR`, `GA`, `Ggtr`, `Gles`.
"""
@relation :keldysh KeldyshCausality(
    GR::RetardedGreensFunction,
    GA::AdvancedGreensFunction,
    Ggtr::GreaterGreensFunction,
    Gles::LesserGreensFunction,
) = (GR - GA) - (Ggtr - Gles)

"""
    AdvancedRetardedConjugate <: AbstractRelation

The advanced propagator is the adjoint of the retarded one,
`G^A(ω) = conj(G^R(ω))` (scalar; `(G^R)†` in orbital space).  So
`G^R − G^A = 2i Im G^R` and the spectral weight is real.

Variables: `GA`, `GR`.  (Complex-valued residual.)
"""
@relation :keldysh AdvancedRetardedConjugate(
    GA::AdvancedGreensFunction, GR::RetardedGreensFunction
) = GA - conj(GR)

# `ζ` names the exchange statistics, so it takes two values and no others. Left
# unconstrained it absorbs whatever ratio the data happens to have: setting it to
# `G^</(e^{-βω} G^>)` zeroes the residual for any pair, so the KMS condition can
# never fail and a pass reports nothing. Measured before this guard: ζ = 7.3 passed.
function _require_statistics_sign(what::Symbol, ζ)
    ζ == 1 ||
        ζ == -1 ||
        error(
            "$what: ζ is the exchange-statistics sign, +1 (bosonic) or -1 (fermionic); " *
            "got $ζ. A free ζ makes the residual vanish for any G^≷ pair, so the " *
            "condition would pass without being tested.",
        )
    return nothing
end

"""
    KeldyshFDT <: AbstractRelation

The fluctuation–dissipation theorem in Keldysh form: in equilibrium the
Keldysh component is fixed by the spectral part through the distribution
function `h`,

`G^K(ω) = h(ω) · (G^R(ω) − G^A(ω))`,

with `h = coth(βω/2)` (bosons) or `tanh(βω/2)` (fermions), supplied by
[`keldysh_distribution`](@ref).  Fluctuation (`G^K`) on the left,
dissipation (`G^R − G^A ∝ Im G^R`) on the right — the two are not
independent in thermal equilibrium.

Supplied-`h` convention, and the caveat that comes with it: `h` here is whatever
the caller says the distribution is, so this tests the FDT only when `h` was
obtained independently. Setting `h = G^K/(G^R − G^A)` zeroes the residual for
ANY state, a fully non-equilibrium one included, since that is the definition of
`h` rather than a law about it. [`KeldyshDistributionEquilibrium`](@ref) is the
half this one does not carry: it pins `h` to `keldysh_distribution`, and the two
together are the theorem.

Variables: `GK`, `h`, `GR`, `GA`.
"""
@relation :keldysh KeldyshFDT(
    GK::KeldyshGreensFunction, h, GR::RetardedGreensFunction, GA::AdvancedGreensFunction
) = GK - h * (GR - GA)

"""
    KeldyshDistributionEquilibrium <: AbstractRelation

The equilibrium value of the Keldysh distribution function,

`h(ω) = coth(βω/2)` (bosonic) or `tanh(βω/2)` (fermionic),

i.e. [`keldysh_distribution`](@ref) written as a relation so a supplied `h` can
be checked rather than assumed.

[`KeldyshFDT`](@ref) alone cannot do this. It reads `h` as given, so it holds out
of equilibrium by construction; this one is what makes "the system is thermal" a
falsifiable claim about the same `h`.

Variables: `h`, `ω`, `β` (or `T`), `stat` ([`Fermionic`](@ref) / [`Bosonic`](@ref)).
"""
@relation :keldysh KeldyshDistributionEquilibrium(h, ω, β::InverseTemperature, stat) =
    h - keldysh_distribution(stat, ω; β=β)

"""
    KMSGreaterLesser <: AbstractRelation

The Kubo–Martin–Schwinger / detailed-balance relation between the
equilibrium greater and lesser correlators,

`G^<(ω) = ζ e^{−βω} G^>(ω)`,   `ζ = +1` (bosons), `ζ = −1` (fermions).

This is the root of the Keldysh FDT: with the RAK identities it forces
`G^K/(G^R − G^A) = (1 + ζe^{−βω})/(1 − ζe^{−βω})`, which is exactly
`coth(βω/2)` (`ζ=+1`) or `tanh(βω/2)` (`ζ=−1`) — see
[`keldysh_distribution`](@ref) and [`KeldyshFDT`](@ref).  Mirrors the
structure-factor [`DetailedBalance`](@ref) `S(−ω) = e^{−βω} S(ω)` on the
propagator side.

Variables: `Gles`, `Ggtr`, `ζ`, `ω`, and `β` (or `T`).
"""
@relation :keldysh KMSGreaterLesser(
    Gles::LesserGreensFunction, Ggtr::GreaterGreensFunction, ζ, ω, β::InverseTemperature
) = begin
    _require_statistics_sign(:KMSGreaterLesser, ζ)
    Gles - ζ * exp(-β * ω) * Ggtr
end

"""
    SpectralFromKeldysh <: AbstractRelation

The bridge between the Keldysh RAK components and the normalized spectral
function `A = −Im G^R/π` (`∫A dω = 1`):

`A(ω) = i (G^R(ω) − G^A(ω)) / (2π)`.

Reduces to [`SpectralFromGreens`](@ref) `A = −Im G^R/π` once
`G^A = conj(G^R)` ([`AdvancedRetardedConjugate`](@ref)) is imposed, so the
real-time and Matsubara spectral definitions agree.

Variables: `A`, `GR`, `GA`.  (Complex-valued residual off equilibrium.)
"""
@relation :keldysh SpectralFromKeldysh(
    A::SpectralFunction, GR::RetardedGreensFunction, GA::AdvancedGreensFunction
) = A - im * (GR - GA) / (2 * π)

# ─── Non-equilibrium: the Langreth rules + the Keldysh kinetic (Dyson) form ───
#
# Out of equilibrium the RAK components of a CONTOUR PRODUCT `C = A·B` follow the
# Langreth rules — exact algebraic identities among the components, equilibrium or
# not, matrix-valued in orbital space (written with `*` so ONE identity holds for
# scalar and matrix propagators alike, like [`Dyson`](@ref)):
#
#   C^R = A^R B^R,   C^A = A^A B^A,
#   C^< = A^R B^< + A^< B^A,   C^> = A^R B^> + A^> B^A.
#
# They preserve the contour structure: a `C` built by these rules satisfies
# [`KeldyshCausality`](@ref) (`C^R − C^A = C^> − C^<`) whenever `A`, `B` do — the
# independent-expectation check in the test.  The steady-state Keldysh Dyson equation
# then fixes the lesser propagator from the lesser self-energy, `G^< = G^R Σ^< G^A`
# (the kinetic form; the equilibrium [`KeldyshFDT`](@ref)/[`KMSGreaterLesser`](@ref)
# are its thermal special case).  The Wigner-transform / gradient-expansion route to
# the Boltzmann equation (Kadanoff–Baym) is a transform → the functional sibling's
# job (scope line, #14).  Reference: Langreth (1976).

"""
    LangrethProductRetarded <: AbstractRelation

The retarded component of a contour product `C = A·B`: `C^R = A^R B^R` (Langreth).
Matrix-valued in orbital space.  Variables: `Cret`, `Aret`, `Bret`.
"""
@relation :keldysh LangrethProductRetarded(Cret, Aret, Bret) = Cret - Aret * Bret

"""
    LangrethProductAdvanced <: AbstractRelation

The advanced component of a contour product `C = A·B`: `C^A = A^A B^A` (Langreth).
Variables: `Cadv`, `Aadv`, `Badv`.
"""
@relation :keldysh LangrethProductAdvanced(Cadv, Aadv, Badv) = Cadv - Aadv * Badv

"""
    LangrethProductLesser <: AbstractRelation

The lesser component of a contour product `C = A·B`, `C^< = A^R B^< + A^< B^A`
(Langreth) — the non-equilibrium generation term.  Matrix-valued.  Variables:
`Cless`, `Aret`, `Bless`, `Aless`, `Badv`.
"""
@relation :keldysh LangrethProductLesser(Cless, Aret, Bless, Aless, Badv) =
    Cless - (Aret * Bless + Aless * Badv)

"""
    LangrethProductGreater <: AbstractRelation

The greater component of a contour product `C = A·B`, `C^> = A^R B^> + A^> B^A`
(Langreth).  Variables: `Cgtr`, `Aret`, `Bgtr`, `Agtr`, `Badv`.
"""
@relation :keldysh LangrethProductGreater(Cgtr, Aret, Bgtr, Agtr, Badv) =
    Cgtr - (Aret * Bgtr + Agtr * Badv)

"""
    KeldyshKineticLesser <: AbstractRelation

The steady-state Keldysh Dyson (kinetic) equation for the lesser propagator,
`G^<(ω) = G^R(ω) Σ^<(ω) G^A(ω)` — the lesser self-energy `Σ^<` drives `G^<` (the
occupation), the non-equilibrium generalization of the FDT lock; in equilibrium it
reduces to [`KeldyshFDT`](@ref)/[`KMSGreaterLesser`](@ref).  Matrix-valued (a triple
product).  `Sless = Σ^<` is a supplied component (not a distinct named quantity).

Variables: `Gless` (`G^<`), `GR` (`G^R`), `Sless` (`Σ^<`), `GA` (`G^A`).
"""
@relation :keldysh KeldyshKineticLesser(
    Gless::LesserGreensFunction,
    GR::RetardedGreensFunction,
    Sless,
    GA::AdvancedGreensFunction,
) = Gless - GR * Sless * GA

"""
    KeldyshKineticGreater <: AbstractRelation

The steady-state Keldysh Dyson (kinetic) equation for the greater propagator,
`G^>(ω) = G^R(ω) Σ^>(ω) G^A(ω)` — the greater partner of
[`KeldyshKineticLesser`](@ref); together `G^≷` fix the occupation and the
in/out-scattering rates.  Matrix-valued.  `Sgtr = Σ^>` is a supplied component.

Variables: `Ggtr` (`G^>`), `GR` (`G^R`), `Sgtr` (`Σ^>`), `GA` (`G^A`).
"""
@relation :keldysh KeldyshKineticGreater(
    Ggtr::GreaterGreensFunction,
    GR::RetardedGreensFunction,
    Sgtr,
    GA::AdvancedGreensFunction,
) = Ggtr - GR * Sgtr * GA

# ─── the self-energy RAK triple + the non-equilibrium distribution ───
#
# The self-energy carries the same RAK structure as the propagator; its
# equilibrium fluctuation–dissipation tie `Σ^K = h(ω)(Σ^R − Σ^A)` mirrors
# [`KeldyshFDT`](@ref) (`Σ^R` is the self-energy of the retarded [`Dyson`](@ref)).
# Off equilibrium the Keldysh propagator is parametrized by a DISTRIBUTION matrix
# `F` (generally non-thermal), `G^K = G^R F − F G^A`, reducing to `KeldyshFDT`
# when `F = h(ω)·I`; for a system between two reservoirs the steady-state `F`
# is the broadening-weighted average of the bath occupations.  Reference:
# Rammer & Smith, [RammerSmith1986](@cite); the two-terminal steady state
# is Haug & Jauho, *Quantum Kinetics in Transport and Optics of Semiconductors*.

"""
    SelfEnergyKeldyshFDT <: AbstractRelation

The self-energy fluctuation–dissipation tie: in equilibrium the Keldysh
self-energy is fixed by the retarded/advanced pair through the distribution
function `h`,

`Σ^K(ω) = h(ω) (Σ^R(ω) − Σ^A(ω))`,

with `h = coth(βω/2)` (bosons) or `tanh(βω/2)` (fermions), supplied by
[`keldysh_distribution`](@ref).  The self-energy mirror of [`KeldyshFDT`](@ref)
(Rammer & Smith, [RammerSmith1986](@cite)).

Variables: `SigmaK` (`Σ^K`), `h`, `SigmaR` (`Σ^R`), `SigmaA` (`Σ^A`).
"""
@relation :keldysh SelfEnergyKeldyshFDT(
    SigmaK::KeldyshSelfEnergy, h, SigmaR::RetardedSelfEnergy, SigmaA::AdvancedSelfEnergy
) = SigmaK - h * (SigmaR - SigmaA)

"""
    NonequilibriumDistribution <: AbstractRelation

The non-equilibrium parametrization of the Keldysh propagator by a DISTRIBUTION
matrix `F` (the generalized, generally non-thermal, occupation),

`G^K(ω) = G^R(ω) F(ω) − F(ω) G^A(ω)`.

This separates the spectral content (`G^R`, `G^A`) from the occupation (`F`); it
reduces to the equilibrium [`KeldyshFDT`](@ref) `G^K = h(G^R − G^A)` when
`F = h(ω)·I` is the thermal scalar distribution.  Matrix-valued in orbital space
(the ordering of the products is kept — `F` need not commute with `G^{R,A}`).
`F = Fdist` is a supplied component.  (Rammer & Smith, [RammerSmith1986](@cite).)

Variables: `GK` (`G^K`), `GR` (`G^R`), `Fdist` (`F`), `GA` (`G^A`).
"""
@relation :keldysh NonequilibriumDistribution(
    GK::KeldyshGreensFunction, GR::RetardedGreensFunction, Fdist, GA::AdvancedGreensFunction
) = GK - (GR * Fdist - Fdist * GA)

"""
    TwoTerminalDistribution <: AbstractRelation

The steady-state non-equilibrium distribution of a region coupled to two
reservoirs `L`, `R` with level-broadenings `Γ_L`, `Γ_R` and occupations `f_L`,
`f_R`,

`f(ω) = [Γ_L(ω) f_L(ω) + Γ_R(ω) f_R(ω)] / [Γ_L(ω) + Γ_R(ω)]`,

the broadening-weighted average of the bath distributions — the non-thermal
occupation that drives a current when `f_L ≠ f_R` (a bias), and the concrete
`F` of [`NonequilibriumDistribution`](@ref).  Reduces to the common equilibrium
occupation when `f_L = f_R`.  (Haug & Jauho, *Quantum Kinetics in Transport and
Optics of Semiconductors*.)

Variables: `fdist` (`f`), `ΓL`, `fL`, `ΓR`, `fR`.
"""
@relation :keldysh TwoTerminalDistribution(fdist, ΓL, fL, ΓR, fR) =
    fdist - (ΓL * fL + ΓR * fR) / (ΓL + ΓR)
