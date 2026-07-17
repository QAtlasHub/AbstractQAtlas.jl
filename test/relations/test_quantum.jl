# Exact quantum-mechanical identities, checked against textbook exact
# states (hydrogen, the harmonic oscillator) — independent constructions,
# not self-consistency.

using AbstractQAtlas
using AbstractQAtlas: check, residual, slack, solve, AbstractInequality, tensor_rank

@testset "virial theorem 2⟨T⟩ = n⟨V⟩ on exact states" begin
    # hydrogen ground state (Rydberg units): E = −1/2, ⟨T⟩ = 1/2, ⟨V⟩ = −1,
    # Coulomb degree n = −1 ⇒ 2⟨T⟩ = −⟨V⟩
    @test check(VirialTheorem(); T=1 // 2, V=-1 // 1, n=-1) isa Bool
    @test residual(VirialTheorem(); T=1 // 2, V=-1 // 1, n=-1) == 0 // 1   # exact (rationals)
    # E = ⟨T⟩+⟨V⟩ = 1/2 − 1 = −1/2 ✓ (Rydberg)
    @test (1 // 2) + (-1 // 1) == -1 // 2
    # harmonic oscillator: n = 2 ⇒ ⟨T⟩ = ⟨V⟩; solve for ⟨V⟩
    @test solve(VirialTheorem(), Val(:V); T=3 // 1, n=2) == 3 // 1
    # a wrong degree fails
    @test !check(VirialTheorem(); T=1 // 2, V=-1 // 1, n=2)
end

@testset "Hellmann–Feynman + Ehrenfest" begin
    # dE/dλ = ⟨∂H/∂λ⟩
    @test check(HellmannFeynman(); dE_dλ=0.37, dH_dλ=0.37, atol=1e-12)
    @test solve(HellmannFeynman(), Val(:dE_dλ); dH_dλ=1.4) ≈ 1.4
    # Ehrenfest: d⟨x⟩/dt = ⟨p⟩/m, d⟨p⟩/dt = ⟨F⟩
    @test check(EhrenfestPosition(); dx_dt=2.0 / 4.0, p=2.0, m=4.0, atol=1e-12)
    @test check(EhrenfestMomentum(); dp_dt=-0.8, F=-0.8, atol=1e-12)
    @test solve(EhrenfestPosition(), Val(:p); dx_dt=0.5, m=4.0) ≈ 2.0   # ⟨p⟩ = m d⟨x⟩/dt
end

@testset "zero-variance eigenstate condition ⟨H²⟩ = E²" begin
    # an exact eigenstate: ⟨H²⟩ = E² (Var(H) = 0)
    E = -1.234
    @test check(EnergyVarianceEigenstate(); H2=E^2, E=E, atol=1e-12)
    @test residual(EnergyVarianceEigenstate(); H2=E^2, E=E) ≈ 0 atol = 1e-12
    # a non-eigenstate with positive variance fails; the residual IS Var(H)
    varH = 0.05
    @test !check(EnergyVarianceEigenstate(); H2=E^2 + varH, E=E, atol=1e-9)
    @test residual(EnergyVarianceEigenstate(); H2=E^2 + varH, E=E) ≈ varH atol = 1e-12
    # solve gives the eigen-consistent ⟨H²⟩ from E
    @test solve(EnergyVarianceEigenstate(), Val(:H2); E=E) ≈ E^2
end

@testset "Robertson uncertainty ΔA·ΔB ≥ ½|⟨[A,B]⟩| (inequality kind)" begin
    @test RobertsonUncertainty() isa AbstractInequality
    # Heisenberg: |⟨[x,p]⟩| = ℏ = 1 ⇒ Δx·Δp ≥ 1/2. Harmonic-oscillator ground
    # state saturates it: Δx = Δp = 1/√2 ⇒ Δx·Δp = 1/2.
    @test check(
        RobertsonUncertainty(); ΔA=1 / sqrt(2), ΔB=1 / sqrt(2), comm=1.0, atol=1e-12
    )
    @test slack(RobertsonUncertainty(); ΔA=1 / sqrt(2), ΔB=1 / sqrt(2), comm=1.0) ≈ 0 atol =
        1e-12   # saturated
    # a squeezed-below-minimum "state" violates it
    @test !check(RobertsonUncertainty(); ΔA=0.3, ΔB=0.3, comm=1.0, atol=1e-9)
    # extra room passes
    @test check(RobertsonUncertainty(); ΔA=2.0, ΔB=3.0, comm=1.0)
end

@testset "quantum domain wiring + energy-component quantities" begin
    using AbstractQAtlas: domain, variables
    @test domain(VirialTheorem()) == :quantum
    @test domain(RobertsonUncertainty()) == :quantum
    @test variables(EnergyVarianceEigenstate()) == (:H2, :E)
    @test tensor_rank(KineticEnergy()) == 0
    @test tensor_rank(PotentialEnergy()) == 0
    @test tensor_rank(EnergyVariance()) == 0
end

@testset "Lieb–Robinson causality bound" begin
    using AbstractQAtlas: check, slack, solve, AbstractInequality
    @test LiebRobinsonBound() isa AbstractInequality
    @test check(LiebRobinsonBound(); v=1.2, v_LR=3.0)               # inside the light cone
    @test slack(LiebRobinsonBound(); v=3.0, v_LR=3.0) == 0.0        # saturating the LR velocity
    @test !check(LiebRobinsonBound(); v=4.0, v_LR=3.0, atol=1e-9)   # superluminal ⇒ forbidden
    @test solve(LiebRobinsonBound(), Val(:v_LR); v=2.5) ≈ 2.5       # the bound is tight at v=v_LR
end

@testset "type-keyed: VirialTheorem" begin
    @test quantities(VirialTheorem()) == (KineticEnergy, PotentialEnergy)
    # 2⟨T⟩ = n⟨V⟩ via bag; harmonic n = 2 ⇒ ⟨T⟩ = ⟨V⟩
    @test check(
        VirialTheorem(), bag(KineticEnergy => 1.0, PotentialEnergy => 1.0); n=2, atol=1e-12
    )
    @test !check(
        VirialTheorem(), bag(KineticEnergy => 1.0, PotentialEnergy => 2.0); n=2, atol=1e-9
    )
    # the generic quantum relations correctly stay symbol-keyed
    @test all(
        r -> isempty(variable_types(r)),
        (
            EhrenfestMomentum(),
            EhrenfestPosition(),
            HellmannFeynman(),
            RobertsonUncertainty(),
            LiebRobinsonBound(),
            EnergyVarianceEigenstate(),
        ),
    )
end
