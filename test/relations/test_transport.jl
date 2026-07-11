# The linear-transport relations, checked against independent textbook
# values and cross-relation consistency.

using AbstractQAtlas
using AbstractQAtlas: check, solve, residual, variables, domain

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

@testset "domain + variables wiring" begin
    @test domain(WiedemannFranz()) == :transport
    @test domain(OpticalSumRule()) == :transport
    @test variables(WiedemannFranz()) == (:κ, :σ, :T, :L0)
    @test variables(MottFormula()) == (:S, :dlnσ_dε, :T)
    @test variables(KelvinRelation()) == (:Π, :S, :T)
end
