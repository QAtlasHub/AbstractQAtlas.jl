# Harris / Luck / Weinrib-Halperin: does a perturbation change the fixed point?
#
# Three published criteria, each a sign. What is checked here is not that they
# were transcribed — that is one line each — but that they AGREE where they must:
# Luck at the random wandering exponent is Harris in one dimension.

using AbstractQAtlas
using AbstractQAtlas: margin

@testset "relevance :: Harris" begin
    # The clean transverse-field Ising chain, ν₀ = 1 at d = 1 (Iglói–Monthus,
    # below Eq. (10.9)). Disorder is relevant — which is why the random chain
    # does not stay at the Ising fixed point.
    @test relevance(HarrisCriterion(); ν₀=1, d=1) === :relevant
    @test margin(HarrisCriterion(); ν₀=1, d=1) == -1

    # ν₀ = 2 at d = 1 is the borderline the review calls out, and it is exact.
    @test relevance(HarrisCriterion(); ν₀=2 // 1, d=1) === :marginal
    @test relevance(HarrisCriterion(); ν₀=3 // 1, d=1) === :irrelevant

    # Raising d lowers the bar, so a fixed point that fails in 1D can pass in 3D.
    @test relevance(HarrisCriterion(); ν₀=1, d=3) === :irrelevant
    # ...and `atol` is the only thing that widens `:marginal` — a float that is
    # merely close is not marginal by default.
    @test relevance(HarrisCriterion(); ν₀=2 + 1e-9, d=1) === :irrelevant
    @test relevance(HarrisCriterion(); ν₀=2 + 1e-9, d=1, atol=1e-6) === :marginal
end

@testset "relevance :: Luck reduces to Harris for random disorder" begin
    # ω = 1/2 is a random sequence, and 1/(1−1/2) = 2 = 2/d at d = 1. The two
    # criteria are then the same statement, so they must agree for EVERY ν₀ —
    # not at one point.
    for ν₀ in (1 // 2, 1 // 1, 2 // 1, 3 // 1, 7 // 2)
        @test margin(LuckCriterion(); ν₀=ν₀, ω=1 // 2) ==
            margin(HarrisCriterion(); ν₀=ν₀, d=1)
        @test relevance(LuckCriterion(); ν₀=ν₀, ω=1 // 2) ===
            relevance(HarrisCriterion(); ν₀=ν₀, d=1)
    end
    # A fixture that CAN disagree: away from ω = 1/2 they must not coincide.
    @test margin(LuckCriterion(); ν₀=1 // 1, ω=0 // 1) !=
        margin(HarrisCriterion(); ν₀=1 // 1, d=1)

    # Fibonacci has ω = −1: bounded fluctuations, so the bar drops to ν₀ > 1/2
    # and the clean Ising chain's ν₀ = 1 clears it.
    @test relevance(LuckCriterion(); ν₀=1, ω=-1) === :irrelevant
    @test relevance(HarrisCriterion(); ν₀=1, d=1) === :relevant   # ...unlike random
end

@testset "relevance :: Weinrib-Halperin takes the DISORDERED ν" begin
    # ρ > 2/ν leaves the uncorrelated class intact. Note the ν here is the
    # uncorrelated disordered fixed point's, not the clean one Harris reads —
    # the two criteria ask different questions of different exponents.
    @test relevance(WeinribHalperinCriterion(); ν=1, ρ=3) === :irrelevant
    @test relevance(WeinribHalperinCriterion(); ν=1, ρ=1) === :relevant
    @test relevance(WeinribHalperinCriterion(); ν=1 // 1, ρ=2 // 1) === :marginal
    # Slower decay is more correlated, hence more relevant — monotone in ρ.
    ms = [margin(WeinribHalperinCriterion(); ν=2, ρ=r) for r in (0.5, 1.0, 2.0, 4.0)]
    @test issorted(ms) && ms[1] < 0 < ms[end]
end

@testset "relevance :: criteria stay out of the relation registry" begin
    # Their residual's SIGN is the answer, so a sweep that reads a nonzero
    # residual as a failure would report every non-marginal system as broken.
    @test !any(r -> r isa RelevanceCriterion, all_relations())
    @test HarrisCriterion() isa RelevanceCriterion
    @test !(HarrisCriterion() isa AbstractRelation)
end
