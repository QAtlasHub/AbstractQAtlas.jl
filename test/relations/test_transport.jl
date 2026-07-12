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

@testset "carrier transport: mobility / Drude / Einstein / Hall" begin
    e = 1.0                                   # natural units
    n, μ = 3.0, 0.4
    # σ = n e μ
    σ = n * e * μ
    @test check(MobilityConductivity(); σ=σ, n=n, e=e, μ=μ, atol=1e-12)
    @test solve(MobilityConductivity(), Val(:μ); σ=σ, n=n, e=e) ≈ μ

    # Drude mobility μ = e τ / m
    τ, m = 0.8, 2.0
    μ_drude = e * τ / m
    @test check(DrudeMobility(); μ=μ_drude, e=e, τ=τ, m=m, atol=1e-12)

    # CROSS-RELATION: MobilityConductivity ∘ DrudeMobility ⇒ σ = n e² τ / m
    σ_drude = solve(MobilityConductivity(), Val(:σ); n=n, e=e, μ=μ_drude)
    @test σ_drude ≈ n * e^2 * τ / m atol = 1e-12

    # Einstein relation μ = e D β, with the β↔T convention
    D, β = 0.6, 1.3
    μ_ein = e * D * β
    @test check(EinsteinRelation(); μ=μ_ein, e=e, D=D, β=β, atol=1e-12)
    @test solve(EinsteinRelation(), Val(:μ); e=e, D=D, T=1 / β) ≈ μ_ein
    # CROSS-RELATION: Einstein + Drude ⇒ D = k_B T τ / m  (= τ/(mβ))
    D_from = solve(EinsteinRelation(), Val(:D); μ=μ_drude, e=e, β=β)
    @test D_from ≈ τ / (m * β) atol = 1e-12

    # single-band Hall coefficient R_H = 1/(n e)
    @test check(SingleBandHall(); R_H=1 / (n * e), n=n, e=e, atol=1e-12)
    @test solve(SingleBandHall(), Val(:R_H); n=n, e=e) ≈ 1 / (n * e)
    @test !check(SingleBandHall(); R_H=2 / (n * e), n=n, e=e, atol=1e-9)

    # Hall angle tan θ_H = σ_xy/σ_xx
    @test check(HallAngle(); tanθ_H=0.3 / 2.0, σxy=0.3, σxx=2.0, atol=1e-12)
    @test solve(HallAngle(), Val(:tanθ_H); σxy=0.3, σxx=2.0) ≈ 0.15
end

@testset "carrier-transport wiring" begin
    @test domain(MobilityConductivity()) == :transport
    @test domain(SingleBandHall()) == :transport
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
