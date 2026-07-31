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

@testset "quantum speed limits: Mandelstam–Tamm + Margolus–Levitin" begin
    using AbstractQAtlas: check, slack, solve, AbstractInequality
    @test MandelstamTammBound() isa AbstractInequality
    @test MargolusLevitinBound() isa AbstractInequality
    # a two-level equal superposition (|0⟩+|1⟩)/√2 of a gap-ω Hamiltonian reaches an
    # ORTHOGONAL state at τ⊥ = π/ω, where ΔE = ω/2 and E−E₀ = ω/2 — the extremal state
    # that SATURATES both bounds simultaneously (an independent physical target: slack 0)
    ω = 1.0
    τ, ΔE, E_above = π / ω, ω / 2, ω / 2
    @test check(MandelstamTammBound(); τ=τ, ΔE=ΔE, atol=1e-12)
    @test slack(MandelstamTammBound(); τ=τ, ΔE=ΔE) ≈ 0 atol = 1e-12
    @test check(MargolusLevitinBound(); τ=τ, E_above=E_above, atol=1e-12)
    @test slack(MargolusLevitinBound(); τ=τ, E_above=E_above) ≈ 0 atol = 1e-12
    # a slower evolution has room to spare: slack = τ − π/(2ΔE) = 2π − π = π
    @test slack(MandelstamTammBound(); τ=2π, ΔE=0.5) ≈ π
    # an impossibly fast orthogonalization (τ below the bound) is forbidden
    @test !check(MandelstamTammBound(); τ=1.0, ΔE=0.5, atol=1e-9)      # 1 < π/(2·0.5) = π
    @test !check(MargolusLevitinBound(); τ=1.0, E_above=0.5, atol=1e-9)
    # solve returns the saturating (minimal) time τ_min = π/(2ΔE)
    @test solve(MandelstamTammBound(), Val(:τ); ΔE=0.5) ≈ π
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
    # The genuinely GENERIC quantum relations stay symbol-keyed, and should: their
    # variables do not name quantities. `RobertsonUncertainty(ΔA, ΔB, comm)` is about
    # two arbitrary observables; `MandelstamTammBound(τ, ΔE)` about a time and an
    # energy spread; `EhrenfestMomentum`/`Position` about expectation values of
    # whichever operator. Typing those would need quantities that do not exist and
    # arguably should not.
    #
    # `LiebRobinsonBound` used to be in this list and is not any more. Its `v_LR` is
    # not a generic symbol — it IS a named quantity, `LiebRobinsonVelocity` — so
    # typing it makes the inequality discoverable from the quantity
    # (`relations_constraining`), which is the whole point of the type-keyed front
    # door. Its other slot, `v`, stays untyped: that is a measured information
    # velocity, and nothing names it yet.
    @test all(
        r -> isempty(variable_types(r)),
        (
            EhrenfestMomentum(),
            EhrenfestPosition(),
            HellmannFeynman(),
            RobertsonUncertainty(),
            MandelstamTammBound(),
            MargolusLevitinBound(),
            EnergyVarianceEigenstate(),
        ),
    )
end

@testset "LiebRobinsonBound is keyed on LiebRobinsonVelocity" begin
    # `v_LR` is a typed subject, so the inequality is discoverable FROM the quantity —
    # `quantities` is typed-subjects ∪ also_constrains, so no manual link is needed.
    # Before it was typed, the relation constrained no named quantity at all, and an
    # atlas holding a v_LR could not find the statement that bounds it.
    # ...and it is no longer in the symbol-keyed group above, which is asserted there.
    @test !isempty(variable_types(LiebRobinsonBound()))
    @test LiebRobinsonVelocity in quantities(LiebRobinsonBound())
    @test any(r -> r isa LiebRobinsonBound, relations_constraining(LiebRobinsonVelocity))

    # the inequality itself is unchanged: slack v_LR - v, non-negative when it holds
    @test residual(LiebRobinsonBound(); v=1.0, v_LR=2.0) == 1.0
    @test check(LiebRobinsonBound(); v=1.0, v_LR=2.0)
    @test !check(LiebRobinsonBound(); v=3.0, v_LR=2.0)
    # saturation is the equality case, and an inequality holds at equality
    @test check(LiebRobinsonBound(); v=2.0, v_LR=2.0)
    @test residual(LiebRobinsonBound(); v=2.0, v_LR=2.0) == 0.0

    # exact arithmetic survives: Rational in, Rational out
    @test residual(LiebRobinsonBound(); v=1 // 2, v_LR=3 // 2) === 1 // 1

    # and it reads from a type-keyed bag through the typed slot
    b = bag(LiebRobinsonVelocity() => 2.0)
    @test haskey(b, AbstractQAtlas._as_key(LiebRobinsonVelocity()))
    @test LiebRobinsonVelocity <: AbstractVelocity
end

@testset "LoschmidtRate: λ = −log L / N ties the rate function to the echo" begin
    # The pair is adopted WITH this law (AbstractQAtlas.jl#128): a definition with
    # only one of its two sides in the vocabulary cannot be stated, which is the
    # "as need arises" that core/quantities.jl's header asks for.
    N = 8
    for L in (1.0, 0.5, 0.1, 1e-6)          # echo ∈ (0, 1]
        λ = -log(L) / N
        @test check(LoschmidtRate(); λ=λ, L=L, N=N, atol=1e-14)
        @test solve(LoschmidtRate(), Val(:λ); L=L, N=N) ≈ λ
    end
    # L = 1 (no decay) is the only zero of the rate function
    @test solve(LoschmidtRate(), Val(:λ); L=1.0, N=N) == 0.0
    @test !check(LoschmidtRate(); λ=0.0, L=0.5, N=N, atol=1e-12)

    # type-keyed: both sides are now real slots, so the pair is CONSTRAINED where
    # before neither name existed here at all
    b = bag(LoschmidtRateFunction => -log(0.5) / N, LoschmidtAmplitude => 0.5)
    @test check(LoschmidtRate(), b; N=N, atol=1e-14)
    @test LoschmidtRate in
        Set(typeof(r) for r in relations_constraining(LoschmidtAmplitude))
    @test LoschmidtRate in
        Set(typeof(r) for r in relations_constraining(LoschmidtRateFunction))

    # `L` enters through log, so the generic solver must REFUSE it rather than
    # linearise — same guard as ξ = v/Δ in Δ.
    err = try
        solve(LoschmidtRate(), Val(:L); λ=0.1, N=N)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("not affine", err.msg)
end
