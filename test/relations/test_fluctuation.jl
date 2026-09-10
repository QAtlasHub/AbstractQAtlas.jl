# test/relations/test_fluctuation.jl — nonequilibrium work fluctuation theorems
# (Jarzynski equality + its second-law corollary, Crooks theorem).

using AbstractQAtlas
using Test
using AbstractQAtlas:
    check, slack, solve, quantities, domain, AbstractInequality, AbstractRelation

@testset "Jarzynski equality + second law (dissipated work ≥ 0)" begin
    # a two-outcome work distribution W ∈ {0, 2} at p = ½, β = 1: compute ⟨e^{−βW}⟩ and
    # ⟨W⟩ directly, read ΔF off Jarzynski, then the second law must hold with STRICTLY
    # positive dissipated work (this is a genuinely irreversible protocol)
    β = 1.0
    exp_work = 0.5 * exp(-β * 0.0) + 0.5 * exp(-β * 2.0)     # ⟨e^{−βW}⟩
    ΔF = -log(exp_work) / β                                  # from ⟨e^{−βW}⟩ = e^{−βΔF}
    W_avg = 0.5 * 0.0 + 0.5 * 2.0                            # ⟨W⟩ = 1.0
    @test JarzynskiEquality() isa AbstractRelation
    @test check(JarzynskiEquality(); exp_work=exp_work, ΔF=ΔF, β=β, atol=1e-12)
    @test JarzynskiSecondLaw() isa AbstractInequality
    @test check(JarzynskiSecondLaw(); W_avg=W_avg, ΔF=ΔF, atol=1e-12)     # ⟨W⟩ ≥ ΔF
    @test slack(JarzynskiSecondLaw(); W_avg=W_avg, ΔF=ΔF) > 0.4           # W_diss ≈ 0.434 > 0
    # β-or-T: passing T instead of β normalizes identically (β = 1 ⇒ T = 1)
    @test check(JarzynskiEquality(); exp_work=exp_work, ΔF=ΔF, T=1.0, atol=1e-12)
    # solve recovers the ΔF-consistent exponential average (affine in exp_work)
    @test solve(JarzynskiEquality(), Val(:exp_work); ΔF=ΔF, β=β) ≈ exp_work
    # a reversible protocol saturates the second law (⟨W⟩ = ΔF ⇒ slack 0)
    @test slack(JarzynskiSecondLaw(); W_avg=0.5, ΔF=0.5) ≈ 0 atol = 1e-12
    # a fabricated second-law violation (⟨W⟩ < ΔF) is caught
    @test !check(JarzynskiSecondLaw(); W_avg=0.3, ΔF=0.5, atol=1e-9)
end

@testset "Crooks fluctuation theorem P_F/P_R = e^{β(W−ΔF)}" begin
    β, ΔF = 1.0, 0.5
    # at the crossing W = ΔF the forward/reverse ratio is exactly 1
    @test check(CrooksFluctuationTheorem(); ratio=1.0, W=ΔF, ΔF=ΔF, β=β, atol=1e-12)
    # off the crossing, ratio = e^{β(W−ΔF)}: W = 1.5 ⇒ e^{1.0}
    @test check(CrooksFluctuationTheorem(); ratio=exp(1.0), W=1.5, ΔF=ΔF, β=β, atol=1e-12)
    # β-or-T normalization applies to the Crooks ratio too
    @test check(CrooksFluctuationTheorem(); ratio=exp(1.0), W=1.5, ΔF=ΔF, T=1.0, atol=1e-12)
    # a fabricated ratio (wrong sign in the exponent) is caught
    @test !check(CrooksFluctuationTheorem(); ratio=exp(-1.0), W=1.5, ΔF=ΔF, β=β, atol=1e-9)
end

@testset "fluctuation domain wiring; supplied-scalar (no named quantity subject)" begin
    @test domain(JarzynskiEquality()) == :fluctuation
    @test domain(JarzynskiSecondLaw()) == :fluctuation
    @test domain(CrooksFluctuationTheorem()) == :fluctuation
    # the averages/ratio are supplied aggregates, not named quantities ⇒ quantities() = ()
    @test all(
        r -> isempty(quantities(r)),
        (JarzynskiEquality(), JarzynskiSecondLaw(), CrooksFluctuationTheorem()),
    )
end

@testset "typical ≤ average, and the gap is how broad the ensemble is" begin
    # Closed-form oracle: for ln X ~ N(μ, σ²) the typical value is exp(μ) and the
    # average is exp(μ + σ²/2), so the RATIO is exp(σ²/2) exactly — known without
    # sampling, and growing without bound in σ.  A broad-distribution fixed point
    # is the σ → ∞ end of this, which is why the two need separate keys.
    μ = -3.0
    for σ in (0.0, 0.5, 1.0, 2.0, 4.0)
        X_typ, X_avg = exp(μ), exp(μ + σ^2 / 2)
        @test check(TypicalBelowAverage(); X_typ=X_typ, X_avg=X_avg)
        @test X_avg / X_typ ≈ exp(σ^2 / 2)
        @test slack(TypicalBelowAverage(); X_typ=X_typ, X_avg=X_avg) ≈ X_avg - X_typ
    end
    # Degenerate distribution (σ = 0) is the saturating case, and only that one.
    @test slack(TypicalBelowAverage(); X_typ=exp(μ), X_avg=exp(μ)) == 0.0
    @test slack(TypicalBelowAverage(); X_typ=exp(μ), X_avg=exp(μ + 0.5)) > 0

    # Independent construction: both reduced from the SAME explicit samples, so
    # the inequality is Jensen on real data rather than on the closed form above.
    rng = MersenneTwister(4242)
    for _ in 1:20
        x = exp.(randn(rng, 200) .* 1.5 .- 2)          # positive, broadly spread
        typ = exp(sum(log, x) / length(x))
        avg = sum(x) / length(x)
        @test check(TypicalBelowAverage(); X_typ=typ, X_avg=avg)
    end

    # ...and it CAN fail: swapping the two reductions is exactly the bug this
    # catches, and a bound nothing can violate would be worth nothing.
    x = exp.(randn(MersenneTwister(7), 200) .* 1.5)
    typ, avg = exp(sum(log, x) / length(x)), sum(x) / length(x)
    @test !check(TypicalBelowAverage(); X_typ=avg, X_avg=typ)

    # The two quantities are DIFFERENT bag keys, which is the point of adding
    # them — a shared key would let one stand in for the other silently.
    @test Typical{MassGap} !== DisorderAveraged{MassGap}
    @test VariableKey(Typical{MassGap}) != VariableKey(DisorderAveraged{MassGap})
    @test VariableKey(Typical{MassGap}) != VariableKey(MassGap)
end
