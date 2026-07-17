# The inequality-relation kind (AbstractInequality) and the entropy /
# quantum-information inequalities, checked against exact bipartite-state
# entropies and their saturating (tight-bound) configurations.

using AbstractQAtlas
using AbstractQAtlas:
    check, residual, slack, solve, AbstractInequality, relation_report, all_relations

@testset "AbstractInequality: check tests the ≥ 0 direction, not |·| ≤ atol" begin
    @test EntropyNonNegativity() isa AbstractInequality
    # satisfied with room, saturated, and violated
    @test check(EntropyNonNegativity(); S=1.3)              # slack 1.3 ≥ 0
    @test check(EntropyNonNegativity(); S=0.0)              # saturated
    @test !check(EntropyNonNegativity(); S=-0.2)            # violated (an equality-style |·| test would wrongly pass at atol≥0.2)
    @test slack(EntropyNonNegativity(); S=1.3) == 1.3
    # a small negative from float noise is tolerated by atol; a real violation is not
    @test check(EntropyNonNegativity(); S=-1e-12, atol=1e-9)
    @test !check(EntropyNonNegativity(); S=-1e-3, atol=1e-9)
    # solve returns the SATURATION value (where the bound is tight)
    @test solve(EntropyNonNegativity(), Val(:S)) == 0
    @test solve(Subadditivity(), Val(:S_AB); S_A=0.6, S_B=0.9) ≈ 1.5   # max S_AB = S_A+S_B
end

@testset "entropy inequalities on a concrete two-qubit pure state" begin
    # a pure entangled state of A⊗B with Schmidt coefficients p, 1−p:
    # S_A = S_B = H(p), S_AB = 0 (pure global state); dims 2.
    H(p) = -p * log(p) - (1 - p) * log(1 - p)
    p = 0.3
    S_A = H(p)
    S_B = H(p)
    S_AB = 0.0
    @test check(Subadditivity(); S_A=S_A, S_B=S_B, S_AB=S_AB, atol=1e-12)      # S_AB ≤ S_A+S_B
    @test check(ArakiLieb(); S_AB=S_AB, S_A=S_A, S_B=S_B, atol=1e-12)          # S_AB ≥ |S_A−S_B| = 0 (saturated)
    @test slack(ArakiLieb(); S_AB=S_AB, S_A=S_A, S_B=S_B) ≈ 0 atol = 1e-12     # a pure state saturates the triangle bound
    @test check(MaxEntropyBound(); S=S_A, log_d=log(2), atol=1e-12)            # S_A ≤ ln 2
    @test check(EntropyNonNegativity(); S=S_A)

    # a maximally entangled (Bell) state saturates MaxEntropyBound
    @test slack(MaxEntropyBound(); S=log(2), log_d=log(2)) ≈ 0 atol = 1e-12
end

@testset "strong subadditivity holds; a fabricated violation is caught" begin
    # SSA: S_AB + S_BC ≥ S_ABC + S_B. Use a valid classical-like assignment
    # (product structure): S_ABC = S_A+S_B+S_C, S_AB=S_A+S_B, S_BC=S_B+S_C.
    S_A, S_B, S_C = 0.4, 0.7, 0.5
    S_AB, S_BC, S_ABC = S_A + S_B, S_B + S_C, S_A + S_B + S_C
    @test check(
        StrongSubadditivity(); S_AB=S_AB, S_BC=S_BC, S_ABC=S_ABC, S_B=S_B, atol=1e-12
    )
    @test slack(StrongSubadditivity(); S_AB=S_AB, S_BC=S_BC, S_ABC=S_ABC, S_B=S_B) ≈ 0 atol =
        1e-12   # product structure saturates SSA
    # break it: inflate S_ABC beyond the bound ⇒ violated
    @test !check(
        StrongSubadditivity(); S_AB=S_AB, S_BC=S_BC, S_ABC=S_ABC + 0.5, S_B=S_B, atol=1e-9
    )
end

@testset "weak monotonicity holds; a fabricated violation is caught" begin
    # WM: S_AB + S_BC ≥ S_A + S_C — the purification dual of SSA. A product structure
    # S_AB = S_A+S_B, S_BC = S_B+S_C gives slack = 2·S_B ≥ 0 (an independent target, not
    # merely ≥0), saturated only when the shared middle S_B = 0.
    S_A, S_B, S_C = 0.4, 0.7, 0.5
    S_AB, S_BC = S_A + S_B, S_B + S_C
    @test check(WeakMonotonicity(); S_AB=S_AB, S_BC=S_BC, S_A=S_A, S_C=S_C, atol=1e-12)
    @test slack(WeakMonotonicity(); S_AB=S_AB, S_BC=S_BC, S_A=S_A, S_C=S_C) ≈ 2 * S_B
    # break it: shrink both pair-unions below S_A + S_C ⇒ violated
    @test !check(WeakMonotonicity(); S_AB=0.3, S_BC=0.3, S_A=S_A, S_C=S_C, atol=1e-9)
end

@testset "entropy of mixing: concavity + Holevo upper bound" begin
    # mixing two ORTHOGONAL pure states equally: S(ρᵢ)=0 ⇒ S_avg=0; the (½,½) mixture has
    # eigenvalues (½,½) ⇒ S_mix=ln2 and H(p)=ln2. Concavity slack = ln2 (mixing raised the
    # entropy from 0); Holevo slack = 0 (orthogonal support saturates the upper bound).
    @test check(EntropyMixingConcavity(); S_mix=log(2), S_avg=0.0, atol=1e-12)
    @test slack(EntropyMixingConcavity(); S_mix=log(2), S_avg=0.0) ≈ log(2)
    @test check(HolevoMixingBound(); S_avg=0.0, H_weights=log(2), S_mix=log(2), atol=1e-12)
    @test slack(HolevoMixingBound(); S_avg=0.0, H_weights=log(2), S_mix=log(2)) ≈ 0 atol =
        1e-12

    # mixing IDENTICAL states: the mixture is that state ⇒ S_mix = S_avg. Concavity is
    # saturated (slack 0); the Holevo bound is loose, slack = H(p) > 0.
    @test slack(EntropyMixingConcavity(); S_mix=1.3, S_avg=1.3) ≈ 0 atol = 1e-12
    @test slack(HolevoMixingBound(); S_avg=1.3, H_weights=log(2), S_mix=1.3) ≈ log(2)

    # fabricated violations of each bound are caught
    @test !check(EntropyMixingConcavity(); S_mix=0.2, S_avg=0.5, atol=1e-9)   # below the average
    @test !check(
        HolevoMixingBound();
        S_avg=0.5,
        H_weights=log(2),
        S_mix=0.5 + log(2) + 0.3,
        atol=1e-9,
    )   # above average + H(p)
end

@testset "Rényi monotonicity S_α non-increasing in α" begin
    # exact Rényi entropies of Schmidt spectrum {p, 1−p}: S_α decreasing in α
    p = 0.25
    Sα(α) = α == 1 ? -p * log(p) - (1 - p) * log(1 - p) : log(p^α + (1 - p)^α) / (1 - α)
    S0, S1, S2, Sinf = Sα(1e-9), Sα(1.0), Sα(2.0), -log(max(p, 1 - p))
    @test check(RenyiMonotonicity(); S_low=S0, S_high=S1, atol=1e-6)      # S_0 ≥ S_1
    @test check(RenyiMonotonicity(); S_low=S1, S_high=S2, atol=1e-12)     # S_1 ≥ S_2
    @test check(RenyiMonotonicity(); S_low=S2, S_high=Sinf, atol=1e-12)   # S_2 ≥ S_∞
    @test !check(RenyiMonotonicity(); S_low=S2, S_high=S1, atol=1e-9)     # reversed order violates
end

@testset "relation_report scores inequalities by their ≥ 0 direction" begin
    # a data set with entropy values: the applicable inequalities must be
    # judged by direction, not |residual| ≤ atol (a positive slack is a PASS).
    data = (; S=0.9, log_d=log(2))
    rep = relation_report(data; atol=1e-12, domain=:entanglement)
    kinds = Dict(typeof(row.relation) => row for row in rep)
    @test haskey(kinds, EntropyNonNegativity)
    @test haskey(kinds, MaxEntropyBound)
    # EntropyNonNegativity: slack 0.9 > 0, so PASS despite |0.9| ≰ 1e-12
    @test kinds[EntropyNonNegativity].residual ≈ 0.9
    @test kinds[EntropyNonNegativity].pass
    @test kinds[MaxEntropyBound].residual ≈ log(2) - 0.9    # negative slack (S=0.9 > ln2≈0.693)
    @test !kinds[MaxEntropyBound].pass                      # 0.9 > ln 2 ⇒ violated
end
