# The linear-transport relations, checked against independent textbook
# values and cross-relation consistency.

using AbstractQAtlas
using AbstractQAtlas: check, solve, residual, variables, domain, tensor_rank

@testset "Wiedemann–Franz κ = L₀ σ T" begin
    L0 = π^2 / 3                       # Sommerfeld Lorenz number (k_B = e = 1)
    σ, T = 2.5, 1.7
    κ = L0 * σ * T
    @test check(WiedemannFranz(); κ=κ, σ=σ, T=T, L0=L0, atol=1e-12)
    # solve for the Lorenz number recovers the Sommerfeld constant
    @test solve(WiedemannFranz(), Val(:L0); κ=κ, σ=σ, T=T) ≈ L0
    # a non-Sommerfeld ratio fails the Sommerfeld check
    @test !check(WiedemannFranz(); κ=κ, σ=σ, T=T, L0=2.0, atol=1e-6)
end

@testset "Mott thermopower S = −(π²/3) T dlnσ/dε" begin
    T, dlnσ = 0.8, 1.4
    S = -(π^2 / 3) * T * dlnσ
    @test check(MottFormula(); S=S, dlnσ_dε=dlnσ, T=T, atol=1e-12)
    @test solve(MottFormula(), Val(:S); dlnσ_dε=dlnσ, T=T) ≈ S
    # zero energy-dependence ⇒ zero thermopower (particle–hole symmetry)
    @test solve(MottFormula(), Val(:S); dlnσ_dε=0.0, T=T) ≈ 0.0
end

@testset "Kelvin relation Π = T S (Onsager)" begin
    S, T = -0.35, 2.0
    @test check(KelvinRelation(); Π=T * S, S=S, T=T, atol=1e-12)
    @test solve(KelvinRelation(), Val(:Π); S=S, T=T) ≈ T * S
    # CROSS-CHECK: Mott S feeds the Kelvin Π — the two relations compose
    dlnσ, Tc = 1.1, 1.3
    S_mott = solve(MottFormula(), Val(:S); dlnσ_dε=dlnσ, T=Tc)
    Π = solve(KelvinRelation(), Val(:Π); S=S_mott, T=Tc)
    @test Π ≈ -(π^2 / 3) * Tc^2 * dlnσ atol = 1e-12   # Π = T·S = −(π²/3)T² dlnσ/dε
end

@testset "Onsager reciprocity L_μν = L_νμ" begin
    # symmetric (zero field): holds
    @test check(OnsagerReciprocity(); L_μν=0.7, L_νμ=0.7, atol=1e-12)
    # asymmetric (e.g. a Hall/antisymmetric part): fails
    @test !check(OnsagerReciprocity(); L_μν=0.7, L_νμ=-0.7, atol=1e-9)
    @test solve(OnsagerReciprocity(), Val(:L_νμ); L_μν=1.3) ≈ 1.3
end

@testset "optical sum rule: total = πD + regular" begin
    D, Wreg = 0.9, 1.6
    total = π * D + Wreg
    @test check(OpticalSumRule(); sigma_integral=total, D=D, W_reg=Wreg, atol=1e-12)
    # extract the Drude weight from the total and regular parts
    @test solve(OpticalSumRule(), Val(:D); sigma_integral=total, W_reg=Wreg) ≈ D
    # a purely regular (insulating) spectrum has zero Drude weight
    @test solve(OpticalSumRule(), Val(:D); sigma_integral=Wreg, W_reg=Wreg) ≈ 0.0 atol =
        1e-12
end

@testset "Johnson–Nyquist FDT: S^j(ω) = ω coth(βω/2) Re σ(ω)" begin
    β, ω, Reσ = 1.5, 0.8, 2.2
    S_j = ω * coth(β * ω / 2) * Reσ
    @test check(CurrentNoiseFDT(); S_j=S_j, Reσ=Reσ, ω=ω, β=β, atol=1e-12)
    @test solve(CurrentNoiseFDT(), Val(:S_j); Reσ=Reσ, ω=ω, β=β) ≈ S_j
    # the β↔T keyword convention (like the other FDTs)
    @test solve(CurrentNoiseFDT(), Val(:S_j); Reσ=Reσ, ω=ω, T=1 / β) ≈ S_j
    # INDEPENDENT EXPECTATION — classical (Nyquist) limit βω ≪ 1: S^j → 2T Re σ
    βs, ωs = 1e-4, 1e-3               # βω = 1e-7 ≪ 1
    S_small = solve(CurrentNoiseFDT(), Val(:S_j); Reσ=Reσ, ω=ωs, β=βs)
    @test S_small ≈ 2 * (1 / βs) * Reσ rtol = 1e-6      # white Nyquist noise
end

@testset "domain + variables wiring" begin
    @test domain(WiedemannFranz()) == :transport
    @test domain(OpticalSumRule()) == :transport
    @test domain(CurrentNoiseFDT()) == :transport
    @test variables(WiedemannFranz()) == (:κ, :σ, :T, :L0)
    @test variables(MottFormula()) == (:S, :dlnσ_dε, :T)
    @test variables(KelvinRelation()) == (:Π, :S, :T)
end

@testset "carrier transport: mobility definition / Einstein / Hall angle (universal)" begin
    e = 1.0                                   # natural units
    n, μ = 3.0, 0.4
    # σ = n e μ (the DEFINITION of the mobility — model-independent)
    σ = n * e * μ
    @test check(MobilityConductivity(); σ=σ, n=n, e=e, μ=μ, atol=1e-12)
    @test solve(MobilityConductivity(), Val(:μ); σ=σ, n=n, e=e) ≈ μ

    # Einstein relation μ = e D β (universal fluctuation–dissipation), β↔T
    D, β = 0.6, 1.3
    μ_ein = e * D * β
    @test check(EinsteinRelation(); μ=μ_ein, e=e, D=D, β=β, atol=1e-12)
    @test solve(EinsteinRelation(), Val(:μ); e=e, D=D, T=1 / β) ≈ μ_ein

    # Hall angle tan θ_H = σ_xy/σ_xx (universal tensor ratio)
    @test check(HallAngle(); tanθ_H=0.3 / 2.0, σxy=0.3, σxx=2.0, atol=1e-12)
    @test solve(HallAngle(), Val(:tanθ_H); σxy=0.3, σxx=2.0) ≈ 0.15
end

@testset "carrier-transport wiring" begin
    @test domain(MobilityConductivity()) == :transport
    @test domain(EinsteinRelation()) == :transport
    @test variables(EinsteinRelation()) == (:μ, :e, :D, :β)
    # the new scalar quantities are rank-0 observables
    for Q in (
        CarrierDensity,
        Mobility,
        ScatteringTime,
        EffectiveMass,
        DiffusionConstant,
        HallCoefficient,
    )
        @test tensor_rank(Q()) == 0
    end
end

@testset "magnetotransport: σ↔ρ tensor inversion, cyclotron, Righi–Leduc, von Klitzing" begin
    # 2×2 conductivity/resistivity inversion: build ρ from σ, then ρ = σ⁻¹ exactly
    σxx, σxy = 2.0, 1.5
    D = σxx^2 + σxy^2
    ρxx = σxx / D
    ρxy = σxy / D
    @test check(LongitudinalResistivity(); ρxx=ρxx, σxx=σxx, σxy=σxy, atol=1e-12)
    @test check(HallResistivity(); ρxy=ρxy, σxx=σxx, σxy=σxy, atol=1e-12)
    # INDEPENDENT: ρ is genuinely the matrix inverse of σ (antisymmetric Hall tensor)
    σ = [σxx σxy; -σxy σxx]
    ρ = [ρxx -ρxy; ρxy ρxx]                    # inv of [[a,b],[-b,a]] is (1/D)[[a,-b],[b,a]]
    @test ρ * σ ≈ [1.0 0.0; 0.0 1.0] atol = 1e-12

    # dissipationless quantum-Hall limit σ_xx = 0 ⇒ ρ_xy = 1/σ_xy, ρ_xx = 0
    @test solve(HallResistivity(), Val(:ρxy); σxx=0.0, σxy=σxy) ≈ 1 / σxy
    @test solve(LongitudinalResistivity(), Val(:ρxx); σxx=0.0, σxy=σxy) ≈ 0.0

    # cyclotron frequency ω_c = eB/m; ties the Hall angle tan θ = ω_c τ
    e, B, m, τ = 1.0, 0.7, 2.0, 0.9
    ωc = solve(CyclotronFrequency(), Val(:ωc); e=e, B=B, m=m)
    @test ωc ≈ e * B / m
    @test check(HallAngle(); tanθ_H=ωc * τ, σxy=ωc * τ, σxx=1.0, atol=1e-12)   # tan θ = ω_c τ

    # Righi–Leduc: κ_xy = L₀ T σ_xy (Wiedemann–Franz in the Hall channel)
    L0, T = π^2 / 3, 1.4
    @test check(RighiLeduc(); κxy=L0 * T * σxy, L0=L0, T=T, σxy=σxy, atol=1e-12)

    # von Klitzing R_xy = h/(ν e²); natural units h = 2π ⇒ R_xy = 2π/(ν e²)
    h, ν = 2π, 3.0
    @test check(VonKlitzing(); Rxy=h / (ν * e^2), ν=ν, e=e, h=h, atol=1e-12)
    @test solve(VonKlitzing(), Val(:Rxy); ν=ν, e=e, h=h) ≈ h / (ν * e^2)
    @test solve(VonKlitzing(), Val(:ν); Rxy=h / (ν * e^2), e=e, h=h) ≈ ν   # read ν off R_xy

    # the new quantities
    @test tensor_rank(Resistivity(:x, :y)) == 2
    @test tensor_rank(MagneticFluxDensity()) == 0
    @test tensor_rank(FillingFactor()) == 0
end

@testset "thermoelectric figure of merit, power factor, Nernst, Ioffe–Regel" begin
    using AbstractQAtlas: slack, AbstractInequality
    S, σ, T, κ = 0.2, 3.0, 1.5, 0.9
    # ZT = S²σT/κ, PF = S²σ, and ZT = PF·T/κ (cross-consistency)
    PF = S^2 * σ
    ZT = S^2 * σ * T / κ
    @test check(PowerFactor(); PF=PF, S=S, σ=σ, atol=1e-12)
    @test check(ThermoelectricFigureOfMerit(); ZT=ZT, S=S, σ=σ, T=T, κ=κ, atol=1e-12)
    @test ZT ≈ PF * T / κ atol = 1e-12                 # ZT = (power factor)·T/κ
    @test solve(PowerFactor(), Val(:PF); S=S, σ=σ) ≈ PF
    @test solve(ThermoelectricFigureOfMerit(), Val(:ZT); S=S, σ=σ, T=T, κ=κ) ≈ ZT

    # Nernst coefficient ν_N = S_xy / B
    @test check(NernstCoefficient(); νN=0.15 / 0.5, S_xy=0.15, B=0.5, atol=1e-12)
    @test solve(NernstCoefficient(), Val(:S_xy); νN=0.3, B=0.5) ≈ 0.15

    # Mott–Ioffe–Regel inequality k_F ℓ ≥ 1 (coherence criterion)
    @test IoffeRegel() isa AbstractInequality
    @test check(IoffeRegel(); kFℓ=5.0)                 # good metal
    @test slack(IoffeRegel(); kFℓ=1.0) == 0.0          # at the MIR limit (saturated)
    @test !check(IoffeRegel(); kFℓ=0.4, atol=1e-9)     # bad metal — criterion violated
end

@testset "type-keyed transport: concrete component keys" begin
    L0, σ, T = π^2 / 3, 3.0, 1.5
    # Wiedemann–Franz through a bag: κ_xx, σ_xx by concrete component, T by field
    @test check(
        WiedemannFranz(),
        bag(
            ThermalConductivity{(:x, :x)} => L0 * σ * T,
            Conductivity{(:x, :x)} => σ,
            Temperature => T,
        );
        L0=L0,
        atol=1e-12,
    )
    # HallAngle: σ_xy and σ_xx are DISTINCT concrete keys in ONE bag — the parametric
    # payoff (a bare-family key could not tell the two components apart)
    σxx, σxy = 3.0, 0.6
    b = bag(Conductivity{(:x, :y)} => σxy, Conductivity{(:x, :x)} => σxx)
    @test check(HallAngle(), b; tanθ_H=σxy / σxx, atol=1e-12)
    @test !check(HallAngle(), b; tanθ_H=σxx / σxy, atol=1e-12)   # swapped ratio ⇒ fails
    # auto-derived quantities: two Conductivity components de-dup to one family
    @test quantities(HallAngle()) == (Conductivity,)
    @test Set(quantities(WiedemannFranz())) == Set((ThermalConductivity, Conductivity))
    # MottFormula's conductivity enters via the supplied derivative dlnσ/dε, kept in
    # the graph through `also_constrains`; RighiLeduc constrains Conductivity via σxy
    @test Set(quantities(MottFormula())) == Set((Thermopower, Conductivity))
    @test Set(quantities(RighiLeduc())) == Set((ThermalConductivity, Conductivity))
    # β-or-T typed conversion carries into transport (Einstein μ = e D β)
    e, D, β = 1.0, 0.5, 2.0
    @test check(
        EinsteinRelation(),
        bag(Mobility => e * D * β, DiffusionConstant => D, InverseTemperature => β);
        e=e,
        atol=1e-12,
    )
end
