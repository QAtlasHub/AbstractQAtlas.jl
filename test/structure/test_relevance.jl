# What is checked is not that the three criteria were transcribed, which is one
# line each, but that they agree where they must: Luck at the random wandering
# exponent is Harris in one dimension.

using AbstractQAtlas
using AbstractQAtlas: margin
using InteractiveUtils: subtypes

# Declared here, not in src: a criterion that forgot its `margin`.
struct _NoMargin <: RelevanceCriterion end

@testset "relevance :: Harris" begin
    # Clean transverse-field Ising chain, ν₀ = 1 at d = 1 (below Eq. (10.9)).
    # Relevant, which is why the random chain leaves the Ising fixed point.
    @test relevance(HarrisCriterion(); ν₀=1, d=1) === :relevant
    @test margin(HarrisCriterion(); ν₀=1, d=1) == -1

    # ν₀ = 2 at d = 1 is the borderline, exactly.
    @test relevance(HarrisCriterion(); ν₀=2 // 1, d=1) === :marginal
    @test relevance(HarrisCriterion(); ν₀=3 // 1, d=1) === :irrelevant

    @test relevance(HarrisCriterion(); ν₀=1, d=3) === :irrelevant   # the bar drops with d
    # Only `atol` widens `:marginal`; merely close is not marginal.
    @test relevance(HarrisCriterion(); ν₀=2 + 1e-9, d=1) === :irrelevant
    @test relevance(HarrisCriterion(); ν₀=2 + 1e-9, d=1, atol=1e-6) === :marginal
end

@testset "relevance :: Luck reduces to Harris for random disorder" begin
    # ω = 1/2 is random, and 1/(1−1/2) = 2 = 2/d at d = 1, so the two are the same
    # statement and must agree for EVERY ν₀, not at one point.
    for ν₀ in (1 // 2, 1 // 1, 2 // 1, 3 // 1, 7 // 2)
        @test margin(LuckCriterion(); ν₀=ν₀, ω=1 // 2) ==
            margin(HarrisCriterion(); ν₀=ν₀, d=1)
        @test relevance(LuckCriterion(); ν₀=ν₀, ω=1 // 2) ===
            relevance(HarrisCriterion(); ν₀=ν₀, d=1)
    end
    # A fixture that CAN disagree: away from ω = 1/2 they must not coincide.
    @test margin(LuckCriterion(); ν₀=1 // 1, ω=0 // 1) !=
        margin(HarrisCriterion(); ν₀=1 // 1, d=1)

    # Fibonacci, ω = −1: bounded fluctuations drop the bar to ν₀ > 1/2, which the
    # clean chain's ν₀ = 1 clears.
    @test relevance(LuckCriterion(); ν₀=1, ω=-1) === :irrelevant
    @test relevance(HarrisCriterion(); ν₀=1, d=1) === :relevant   # ...unlike random
end

@testset "relevance :: Weinrib-Halperin takes the DISORDERED ν" begin
    # ρ > 2/ν leaves the uncorrelated class intact. This ν is the uncorrelated
    # DISORDERED fixed point's, not the clean ν₀ Harris reads.
    @test relevance(WeinribHalperinCriterion(); ν_dis=1, ρ=3) === :irrelevant
    @test relevance(WeinribHalperinCriterion(); ν_dis=1, ρ=1) === :relevant
    @test relevance(WeinribHalperinCriterion(); ν_dis=1 // 1, ρ=2 // 1) === :marginal
    # Slower decay is more correlated, so more relevant: monotone in ρ.
    ms = [margin(WeinribHalperinCriterion(); ν_dis=2, ρ=r) for r in (0.5, 1.0, 2.0, 4.0)]
    @test issorted(ms) && ms[1] < 0 < ms[end]
end

@testset "relevance :: a non-finite margin is not a verdict" begin
    # IEEE makes both `abs(NaN) <= atol` and `NaN > 0` false, so the else-branch
    # would report `:relevant`. No atol rescues it, not even Inf.
    for kw in ((; ν₀=NaN, d=1.0), (; ν₀=Inf, d=Inf))
        @test_throws ErrorException relevance(HarrisCriterion(); kw...)
        @test_throws ErrorException relevance(HarrisCriterion(); atol=Inf, kw...)
    end

    # The boundaries that produced them are refused at the source.
    @test_throws ArgumentError margin(HarrisCriterion(); ν₀=1, d=0)
    @test_throws ArgumentError margin(HarrisCriterion(); ν₀=1, d=-2)
    # ±0.0 are both "zero dimensions" and gave OPPOSITE verdicts, since 2/(-0.0)
    # is -Inf.
    @test_throws ArgumentError margin(HarrisCriterion(); ν₀=1.0, d=-0.0)
    @test_throws ArgumentError margin(LuckCriterion(); ν₀=1, ω=1)       # divides by zero
    @test_throws ArgumentError margin(LuckCriterion(); ν₀=1, ω=2)       # outgrows L
    @test_throws ArgumentError margin(WeinribHalperinCriterion(); ν_dis=0, ρ=1)
    @test_throws ArgumentError margin(WeinribHalperinCriterion(); ν_dis=1, ρ=0)
end

@testset "relevance :: every criterion answers the interface" begin
    # As test_invariants.jl sweeps subtypes: a criterion added later without a
    # `margin` fails HERE, not as a MethodError somewhere downstream.
    subs = subtypes(RelevanceCriterion)
    @test length(subs) >= 3
    for T in subs
        @test hasmethod(margin, Tuple{T})            # ...and not just the fallback
    end

    # ...and the fallback names what is missing.
    msg = try
        margin(_NoMargin(); ν₀=1, d=1)
        ""
    catch err
        err isa ErrorException ? sprint(showerror, err) : rethrow()
    end
    @test occursin("RelevanceCriterion", msg) && occursin("no `margin` method", msg)
end

@testset "relevance :: criteria stay out of the relation registry" begin
    # The SIGN is the answer, so a sweep reading a nonzero residual as a failure
    # would call every non-marginal system broken.
    @test !any(r -> r isa RelevanceCriterion, all_relations())
    @test HarrisCriterion() isa RelevanceCriterion
    @test !(HarrisCriterion() isa AbstractRelation)
end
