# Entanglement-entropy relations vs INDEPENDENT constructions:
# purity from an explicit density matrix, the cut-counting coefficient read off
# the EXACT free-fermion critical Ising chain, and Page's formula against exact
# small cases and a Haar-random-state average.
#
# The free-fermion arm is here rather than a synthetic `S(ℓ) = (c/3) ln ℓ`
# because a fixture built from the relation cannot test the relation: it
# differentiates the coefficient it was handed and gets it back, so it reads the
# same whether the coefficient is `c/3`, `c/6` or `ncuts·c/6`, and it cannot see
# the geometry at all.  The chain below is quadratic in Majoranas, so `S(ℓ)` is
# exact at any size, and the END block and the CENTRED block are measured on the
# SAME ground state — leaving the cut count as the only difference between them.

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

@testset "CFT entanglement slope counts cuts: exact critical Ising chain" begin
    # Fixture: the CLEAN critical transverse-field Ising chain with OPEN ends,
    # H = -Σ Z_j Z_{j+1} - Σ X_j, which is c = 1/2 and quadratic in Majoranas
    # a_{2j-1}, a_{2j}:  H = (i/4) Σ A_{mn} a_m a_n  with  A[2j-1,2j] = -2h_j,
    # A[2j,2j+1] = -2J_j.  The ground-state covariance is the orthogonal polar
    # factor of A (each ε_k in the canonical form replaced by 1, orientation
    # kept), and the entropy of a region is Peschel's function of the restricted
    # covariance spectrum — i.e. the package's own
    # `free_fermion_entanglement_entropy`.
    c = 1 / 2

    function covariance(N)
        A = zeros(2N, 2N)
        for j in 1:N
            A[2j - 1, 2j] = -2.0        # -2 h_j, h_j = 1
        end
        for j in 1:(N - 1)
            A[2j, 2j + 1] = -2.0        # -2 J_j, J_j = 1
        end
        A = A - transpose(A)            # antisymmetry as structure, not as bookkeeping
        F = svd(A)
        Γ = F.U * F.Vt
        return (Γ .- transpose(Γ)) ./ 2   # antisymmetry drifts at a zero singular value
    end

    function entropy(Γ, sites)
        idx = vcat(([2j - 1, 2j] for j in sites)...)
        ν = eigvals(Hermitian(im .* Γ[idx, idx]))
        return free_fermion_entanglement_entropy([
            (1 + real(x)) / 2 for x in ν if real(x) > 0
        ])
    end

    # OLS slope of S against ln ℓ; the additive constant drops out of the weights
    function log_slope(ℓs, S)
        x = log.(float.(ℓs))
        x .-= sum(x) / length(x)
        return sum(x .* S) / sum(abs2, x)
    end

    # ONE cut: the block sits at an open end, so only its inner edge is a cut.
    # TWO cuts: the block sits in the bulk, both edges away from the boundary.
    end_block(N, ℓ) = 1:ℓ
    centred_block(N, ℓ) = (div(N - ℓ, 2) + 1):(div(N - ℓ, 2) + ℓ)

    measured = map((48, 128)) do N
        Γ = covariance(N)
        ℓs = 3:div(N, 4)          # short of ℓ ≈ N, where a centred block fills the chain
        a1 = log_slope(ℓs, [entropy(Γ, end_block(N, ℓ)) for ℓ in ℓs])
        a2 = log_slope(ℓs, [entropy(Γ, centred_block(N, ℓ)) for ℓ in ℓs])
        (; N, a1, a2)
    end

    for m in measured
        # The relation, on exact data, at both cut counts.  `check` takes an
        # ABSOLUTE tolerance, and each is the MEASURED deviation rounded up by a
        # factor of about three, so a regression well below the finite-size floor
        # still fails:
        #
        #   ncuts = 1:  0.000166 (N = 48), 0.000179 (N = 128)  — 0.20-0.21 % of c/6
        #   ncuts = 2:  0.011596 (N = 48), 0.007269 (N = 128)  — 6.96 → 4.36 % of c/3
        #
        # The two-cut region carries the larger correction because BOTH of its
        # edges sit a finite distance from the ends of the chain, and it shrinks
        # with N, as it must.  The fixture is deterministic to ~1e-14, so these
        # are not noise budgets.
        @test check(CFTEntanglementSlope(); dS_dlogℓ=m.a1, c=c, ncuts=1, atol=0.0005)
        @test check(CFTEntanglementSlope(); dS_dlogℓ=m.a2, c=c, ncuts=2, atol=0.015)

        # ...and it can DISAGREE: each geometry is nearer its own cut count than
        # the other's.  Without this the two `check`s above would both pass a
        # relation that ignored `ncuts` entirely.
        @test abs(m.a1 - 1 * c / 6) < abs(m.a1 - 2 * c / 6)
        @test abs(m.a2 - 2 * c / 6) < abs(m.a2 - 1 * c / 6)

        # `c` cancels in the ratio, so this is `ncuts` and nothing else.  Measured
        # |ratio − 2| = 0.143 (N = 48), 0.091 (N = 128); the budget is the same
        # ~1.3× round-up as the two-cut line above.  Discriminating power is against
        # ONE cut, which would put the ratio at 1, a whole unit away.
        @test m.a2 / m.a1 ≈ 2 atol = 0.18

        # The same measured slope yields a DIFFERENT central charge under a
        # different cut count — which is why `ncuts` is required, not defaulted.
        @test solve(CFTEntanglementSlope(), Val(:c); dS_dlogℓ=m.a2, ncuts=2) ≈ c rtol = 0.08
        @test solve(CFTEntanglementSlope(), Val(:c); dS_dlogℓ=m.a2, ncuts=1) ≈ 2c rtol =
            0.08
    end

    # Vary an axis that must not matter: the coefficient is a property of the
    # fixed point, so a larger chain moves the ratio TOWARDS 2, never away.
    @test abs(measured[2].a2 / measured[2].a1 - 2) <
        abs(measured[1].a2 / measured[1].a1 - 2)

    # Pure arithmetic on the declared form, with no fixture in the way.
    @test check(CFTEntanglementSlope(); dS_dlogℓ=1 / 3, c=1.0, ncuts=2, atol=1e-14)
    @test check(CFTEntanglementSlope(); dS_dlogℓ=1 / 6, c=1.0, ncuts=1, atol=1e-14)

    # Auto-discovery consequence, and the reason `ncuts` is a declared variable
    # rather than a default: data that does not say which geometry the slope was
    # measured on yields NO central-charge row, instead of one that assumed a
    # geometry on the caller's behalf.
    @test isempty(relation_report((; dS_dlogℓ=1 / 3, c=1.0)))
    discovered = relation_report((; dS_dlogℓ=1 / 3, c=1.0, ncuts=2))
    @test only(discovered).relation isa CFTEntanglementSlope
    @test only(discovered).pass
end

@testset "CFTEntanglementSlope is type-keyed like its cft.jl siblings" begin
    rel = CFTEntanglementSlope()

    # `c` is the typed subject, as in CasimirCentralCharge / CardyDensityOfStates.
    @test variable_types(rel) == (CentralCharge,)
    @test CentralCharge in quantities(rel)
    @test rel in relations_constraining(CentralCharge)

    # ...and the entropy is not lost by typing `c`: it enters through the SUPPLIED
    # derivative, so it has no slot of its own and is declared via `also_constrains`
    # — the hook that exists for exactly this shape.
    @test also_constrains(rel) == (VonNeumannEntropy,)
    @test VonNeumannEntropy in quantities(rel)
    @test rel in relations_constraining(VonNeumannEntropy)

    # `ncuts` names no quantity, so it stays a supplied coordinate like `L` in
    # CasimirCentralCharge; it must not have become a bag-visible subject.
    @test :ncuts in variables(rel)
    @test length(variable_types(rel)) == 1
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

@testset "free-fermion entanglement from the correlation matrix (Peschel)" begin
    using AbstractQAtlas:
        check, solve, free_fermion_entanglement_entropy, free_fermion_renyi_entropy
    using LinearAlgebra: eigvals, Hermitian

    # GENUINE Gaussian case: the bonding orbital (c†₁+c†₂)/√2|0⟩ of a 2-site
    # hopping model.  Correlation matrix C_ij = ⟨c†_i c_j⟩ = 1/2 (all entries);
    # trace out site 2 ⇒ C_A = [1/2], eigenvalue ζ = 1/2 ⇒ S_A = ln 2 (one Bell pair).
    C = [0.5 0.5; 0.5 0.5]
    C_A = C[1:1, 1:1]                              # region A = site 1
    ζ = real(eigvals(Hermitian(C_A)))
    @test ζ ≈ [0.5] atol = 1e-12
    @test free_fermion_entanglement_entropy(ζ) ≈ log(2) atol = 1e-12

    # a filled/empty mode contributes nothing; a maximal mode ln 2
    @test free_fermion_entanglement_entropy([0.0, 1.0, 1.0]) == 0.0
    @test free_fermion_entanglement_entropy([0.5, 0.5, 0.5]) ≈ 3 * log(2) atol = 1e-12

    # Peschel single-particle spectrum ε = ln((1−ζ)/ζ): ζ=1/2 ⇒ ε=0 (max entangled)
    @test check(EntanglementSpectrumCorrelation(); ε=0.0, ζ=0.5, atol=1e-12)
    @test solve(EntanglementSpectrumCorrelation(), Val(:ε); ζ=0.3) ≈ log(0.7 / 0.3)
    # inverting the spectrum recovers the Fermi-Dirac occupation ζ = 1/(e^ε+1)
    ε = 1.1
    ζinv = 1 / (exp(ε) + 1)
    @test check(EntanglementSpectrumCorrelation(); ε=ε, ζ=ζinv, atol=1e-12)

    # Rényi → von Neumann as n → 1 (INDEPENDENT: two different formulas agree)
    ζset = [0.15, 0.5, 0.82, 0.97]
    @test free_fermion_renyi_entropy(ζset, 1.0001) ≈ free_fermion_entanglement_entropy(ζset) atol =
        1e-3
    # Rényi-2 = −Σ ln(ζ²+(1−ζ)²) matches the moment definition per mode
    @test free_fermion_renyi_entropy([0.3], 2) ≈ -log(0.3^2 + 0.7^2) atol = 1e-12
    @test_throws ErrorException free_fermion_renyi_entropy(ζset, 1)   # n=1 is the vN limit
end

@testset "multipartite: monogamy / three-tangle on GHZ and W states" begin
    using AbstractQAtlas: check, solve, slack, AbstractInequality

    # tangle = concurrence²; a Bell pair has C=1 ⇒ τ=1
    @test check(ConcurrenceTangle(); τ=1.0, C=1.0, atol=1e-12)
    @test solve(ConcurrenceTangle(), Val(:τ); C=2 / 3) ≈ 4 / 9   # W-state pair concurrence

    # GHZ = (|000⟩+|111⟩)/√2: pairwise reduced states are separable ⇒ τ_AB=τ_AC=0,
    # while A is maximally entangled with BC ⇒ τ(A:BC)=1, so the three-tangle = 1.
    τ_ABC_ghz, τ_AB_ghz, τ_AC_ghz = 1.0, 0.0, 0.0
    @test check(Monogamy(); τ_ABC=τ_ABC_ghz, τ_AB=τ_AB_ghz, τ_AC=τ_AC_ghz)
    τ3_ghz = solve(
        ThreeTangleDefinition(), Val(:τ3); τ_ABC=τ_ABC_ghz, τ_AB=τ_AB_ghz, τ_AC=τ_AC_ghz
    )
    @test τ3_ghz ≈ 1.0 atol = 1e-12                              # GHZ: genuine tripartite entanglement
    @test τ3_ghz ≈ slack(Monogamy(); τ_ABC=τ_ABC_ghz, τ_AB=τ_AB_ghz, τ_AC=τ_AC_ghz) atol =
        1e-12

    # W = (|001⟩+|010⟩+|100⟩)/√3: τ(A:B)=τ(A:C)=4/9, τ(A:BC)=8/9 ⇒ three-tangle = 0
    # (W saturates monogamy — no residual tripartite tangle).
    τ_ABC_w, τ_AB_w, τ_AC_w = 8 / 9, 4 / 9, 4 / 9
    @test check(Monogamy(); τ_ABC=τ_ABC_w, τ_AB=τ_AB_w, τ_AC=τ_AC_w, atol=1e-12)
    @test slack(Monogamy(); τ_ABC=τ_ABC_w, τ_AB=τ_AB_w, τ_AC=τ_AC_w) ≈ 0 atol = 1e-12
    @test solve(
        ThreeTangleDefinition(), Val(:τ3); τ_ABC=τ_ABC_w, τ_AB=τ_AB_w, τ_AC=τ_AC_w
    ) ≈ 0 atol = 1e-12

    # monogamy is a genuine bound: over-sharing (τ_AB+τ_AC > τ_ABC) is forbidden
    @test Monogamy() isa AbstractInequality
    @test !check(Monogamy(); τ_ABC=0.5, τ_AB=0.4, τ_AC=0.4, atol=1e-9)
end

@testset "tripartite information + Kitaev–Preskill TEE" begin
    using AbstractQAtlas: check, solve
    # I₃ = I(A:B) + I(A:C) − I(A:BC)
    @test check(
        TripartiteInformationDefinition();
        I3=0.3 + 0.5 - 1.1,
        I_AB=0.3,
        I_AC=0.5,
        I_ABC=1.1,
        atol=1e-12,
    )
    @test solve(
        TripartiteInformationDefinition(), Val(:I3); I_AB=0.3, I_AC=0.5, I_ABC=1.1
    ) ≈ -0.3

    # Kitaev–Preskill: the alternating tripartite sum isolates −γ. Build a
    # toric-code-like assignment (all singles a, pairs b, triple c):
    # 3a − 3b + c = −γ ⇒ pick a,b,c so γ = ln 2 (toric code total dimension D=2).
    a, b = 1.7, 2.4
    γ = log(2)
    c = -γ - (3a - 3b)                          # so 3a − 3b + c = −γ
    @test check(
        KitaevPreskillTEE();
        γ=γ,
        S_A=a,
        S_B=a,
        S_C=a,
        S_AB=b,
        S_BC=b,
        S_CA=b,
        S_ABC=c,
        atol=1e-12,
    )
    @test solve(
        KitaevPreskillTEE(), Val(:γ); S_A=a, S_B=a, S_C=a, S_AB=b, S_BC=b, S_CA=b, S_ABC=c
    ) ≈ log(2) atol = 1e-12                      # extracts γ = ln 2

    @test tensor_rank(Concurrence()) == 0
    @test tensor_rank(TopologicalEntanglementEntropy()) == 0
end
