# Entanglement-entropy relations vs INDEPENDENT constructions:
# purity from an explicit density matrix, the central charge read off a
# synthetic CFT log-growth, and Page's formula against exact small cases
# and a Haar-random-state average.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve
using LinearAlgebra, Random

@testset "Rényi-2 from purity: S_2 = −ln Tr ρ²" begin
    # pure state: purity 1, S_2 = 0
    @test solve(RenyiTwoPurity(), Val(:S2); purity=1.0) == 0.0
    # maximally mixed on d levels: purity = 1/d, S_2 = ln d
    for d in (2, 4, 8)
        ρ = Matrix{Float64}(I, d, d) ./ d
        purity = tr(ρ^2)                      # = 1/d, built independently
        @test check(RenyiTwoPurity(); S2=log(d), purity=purity, atol=1e-13)
    end
    # a generic ρ: S_2 vs −ln Tr ρ² from the explicit matrix
    ρ = [0.6 0.1; 0.1 0.4]
    @test check(RenyiTwoPurity(); S2=(-log(tr(ρ^2))), purity=tr(ρ^2), atol=1e-13)
end

@testset "CFT entanglement slope reads off the central charge" begin
    # synthetic S(ℓ) = (c/3) ln ℓ + const with c = 1/2 (Ising); the slope
    # in ln ℓ must return c/3, independent of the (dropped) constant.
    c = 1 / 2
    S(ℓ) = (c / 3) * log(ℓ) + 0.77
    ℓ = 40.0
    h = 1e-4
    dS_dlogℓ = (S(ℓ * exp(h)) - S(ℓ * exp(-h))) / (2h)   # d/d(ln ℓ)
    @test check(CFTEntanglementSlope(); dS_dlogℓ=dS_dlogℓ, c=c, atol=1e-6)
    @test solve(CFTEntanglementSlope(), Val(:c); dS_dlogℓ=dS_dlogℓ) ≈ c atol = 1e-6
    # c = 1 free boson: slope 1/3
    @test check(CFTEntanglementSlope(); dS_dlogℓ=1 / 3, c=1.0, atol=1e-14)
end

@testset "Page average entropy: exact small cases + symmetry" begin
    # two qubits (dA=dB=2): ⟨S⟩ = 1/3 + 1/4 − 1/4 = 1/3, exactly known
    @test page_average_entropy(2, 2) ≈ 1 / 3
    # symmetric in the two dimensions
    @test page_average_entropy(2, 64) == page_average_entropy(64, 2)
    # nearly maximal: ⟨S⟩ → ln(dA) as dB → ∞, with a small positive deficit
    # of order dA/(2 dB) (the harmonic-sum correction is the same order, so
    # this is an asymptotic-form sanity, not a tight match).
    @test page_average_entropy(2, 128) < log(2)                 # below maximal
    @test log(2) - page_average_entropy(2, 128) < 2 / (2 * 128) * 1.5   # deficit ~ dA/(2dB)
    @test page_average_entropy(2, 4096) ≈ log(2) rtol = 1e-3     # deficit vanishes as dB→∞
    @test_throws ErrorException page_average_entropy(0, 4)
end

@testset "Page formula vs a Haar-random-state average" begin
    # dA=2, dB=3: sample random pure states, average the reduced-ρ_A entropy
    dA, dB = 2, 3
    rng = MersenneTwister(2024)
    Ssum = 0.0
    M = 4000
    for _ in 1:M
        ψ = randn(rng, ComplexF64, dA * dB)
        ψ ./= norm(ψ)
        Ψ = reshape(ψ, dA, dB)              # A ⊗ B
        ρA = Ψ * Ψ'                         # reduced density matrix on A
        ev = filter(>(1e-14), real(eigvals(Hermitian(ρA))))
        Ssum += -sum(ev .* log.(ev))
    end
    @test isapprox(Ssum / M, page_average_entropy(dA, dB); rtol=0.03)
end

@testset "entropy zoo: Rényi/Tsallis moments, mutual & conditional, Klein" begin
    using AbstractQAtlas: check, solve, slack, residual, AbstractInequality

    # Schmidt spectrum {p, 1−p}: moments Tr ρ^α = p^α + (1−p)^α
    p = 0.3
    mom(α) = p^α + (1 - p)^α
    SvN = -p * log(p) - (1 - p) * log(1 - p)

    # Rényi from the moment reduces to −ln(purity) at α=2 (= RenyiTwoPurity) …
    S2 = solve(RenyiEntropyMoment(), Val(:Sα); moment=mom(2), α=2)
    @test S2 ≈ -log(mom(2)) atol = 1e-12
    @test check(RenyiTwoPurity(); S2=S2, purity=mom(2), atol=1e-12)      # cross-relation
    # … and → S_vN as α → 1 (l'Hôpital limit, numerically)
    @test solve(RenyiEntropyMoment(), Val(:Sα); moment=mom(1.0001), α=1.0001) ≈ SvN atol =
        1e-3

    # Tsallis from the moment; q → 1 limit is S_vN
    @test check(
        TsallisEntropyMoment(); Sq=(1 - mom(2)) / (2 - 1), moment=mom(2), q=2, atol=1e-12
    )
    @test solve(TsallisEntropyMoment(), Val(:Sq); moment=mom(1.0001), q=1.0001) ≈ SvN atol =
        1e-3

    # mutual information I = S_A+S_B−S_AB ≥ 0, and on a pure state (S_AB=0)
    # with S_A=S_B=H(p): I = 2 H(p); consistency with Subadditivity slack
    S_A = S_B = SvN
    S_AB = 0.0
    I = solve(MutualInformationDefinition(), Val(:I); S_A=S_A, S_B=S_B, S_AB=S_AB)
    @test I ≈ 2 * SvN atol = 1e-12
    @test slack(Subadditivity(); S_A=S_A, S_B=S_B, S_AB=S_AB) ≈ I atol = 1e-12   # I = subadditivity slack

    # conditional entropy of a pure entangled state is NEGATIVE: S(A|B)=S_AB−S_B=−S_B
    S_cond = solve(ConditionalEntropyDefinition(), Val(:S_cond); S_AB=S_AB, S_B=S_B)
    @test S_cond ≈ -S_B atol = 1e-12
    @test S_cond < 0                                                     # entanglement witness

    # Klein's inequality S(ρ‖σ) ≥ 0 (inequality kind), zero iff ρ=σ
    @test RelativeEntropyNonNegativity() isa AbstractInequality
    @test check(RelativeEntropyNonNegativity(); S_rel=0.4)
    @test slack(RelativeEntropyNonNegativity(); S_rel=0.0) == 0.0        # ρ = σ saturates
    @test !check(RelativeEntropyNonNegativity(); S_rel=-1e-3, atol=1e-9)
end

@testset "measurement + Markov entropies on concrete states" begin
    using AbstractQAtlas: check, solve, slack, AbstractInequality
    using LinearAlgebra: eigvals, diagm, Hermitian

    # a 2×2 density matrix with coherences; measure (dephase) in the z basis
    a, c = 0.7, 0.3               # populations
    off = 0.35                    # coherence (|off|² ≤ a·c for positivity: 0.1225 ≤ 0.21 ✓)
    ρ = [a off; off c]
    Δρ = [a 0.0; 0.0 c]           # dephased (diagonal part)
    ent(M) = (λ=filter(>(1e-15), real(eigvals(Hermitian(M)))); -sum(x -> x * log(x), λ))
    S = ent(ρ)
    S_meas = ent(Δρ)
    # relative entropy S(ρ‖Δρ) = Tr ρ(ln ρ − ln Δρ)
    using LinearAlgebra: tr
    lnρ = log(Hermitian(ρ))
    lnΔρ = diagm([log(a), log(c)])
    S_rel = real(tr(ρ * (lnρ - lnΔρ)))

    # measurement does not decrease entropy: S(Δρ) ≥ S(ρ)
    @test check(MeasurementEntropyIncrease(); S_meas=S_meas, S=S)
    @test slack(MeasurementEntropyIncrease(); S_meas=S_meas, S=S) > 0     # strict, coherences present
    # …and the gain equals the relative entropy to the dephased state (exact identity)
    @test check(MeasurementEntropyRelative(); S_meas=S_meas, S=S, S_rel=S_rel, atol=1e-10)
    @test (S_meas - S) ≈ S_rel atol = 1e-10
    # a state already diagonal saturates the inequality (Δρ = ρ)
    @test slack(MeasurementEntropyIncrease(); S_meas=S, S=S) == 0.0

    # MarkovEntropy = conditional mutual information = strong-subadditivity slack
    S_A, S_B, S_C = 0.4, 0.7, 0.5
    S_AB, S_BC, S_ABC = S_A + S_B, S_B + S_C, S_A + S_B + S_C          # product (Markov) state
    I_cmi = solve(
        MarkovEntropyDefinition(), Val(:I_cmi); S_AB=S_AB, S_BC=S_BC, S_ABC=S_ABC, S_B=S_B
    )
    @test I_cmi ≈ slack(StrongSubadditivity(); S_AB=S_AB, S_BC=S_BC, S_ABC=S_ABC, S_B=S_B) atol =
        1e-12
    @test I_cmi ≈ 0 atol = 1e-12                                        # a Markov chain: I(A:C|B)=0

    @test MeasurementEntropyIncrease() isa AbstractInequality
    @test tensor_rank(MeasurementEntropy()) == 0
    @test tensor_rank(MarkovEntropy()) == 0
end
