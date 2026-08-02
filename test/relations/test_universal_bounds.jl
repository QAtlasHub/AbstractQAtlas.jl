# The eight universal bounds: each states `measured OP fetched-bounding-value`,
# with the bounding side typed and the bounded side an untyped slot.
#
# The numbers below are the closed forms the bounding quantities carry, written
# out here rather than fetched — this package holds no `fetch`, so the physics
# check is that the STATEMENT has the right direction and saturates where the
# literature says it does, not that some atlas returns the right number.

using AbstractQAtlas
using Test
using AbstractQAtlas:
    AbstractInequality,
    bound_direction,
    bounded_slot,
    bounding_slot,
    bounds_on,
    check,
    quantities,
    slack,
    solve,
    variable_slots

# (relation, bounded, bounding, direction, bounding quantity)
const _EIGHT = [
    (CHSHInequality(), :S, :S_max, :upper, CHSHBound),
    (MerminInequality(), :M, :M_max, :upper, MerminGHZBound),
    (LyapunovChaosBound(), :λ_L, :λ_max, :upper, ChaosBound),
    (OrthogonalizationTimeBound(), :τ, :τ_min, :lower, QuantumSpeedLimit),
    (FastScramblingBound(), :t_scr, :t_min, :lower, ScramblingTime),
    (SecretKeyRateBound(), :r, :r_min, :lower, BB84KeyRate),
    (CloningFidelityBound(), :F, :F_max, :upper, OptimalCloningFidelity),
    (BekensteinEntropyBound(), :S, :S_max, :upper, BekensteinBound),
]

@testset "each bound types its BOUNDING side and leaves the bounded side open" begin
    for (rel, bd, bg, dir, Q) in _EIGHT
        @test rel isa AbstractInequality
        @test bounded_slot(rel) === bd
        @test bounding_slot(rel) === bg
        @test bound_direction(rel) === dir
        # the asymmetry that makes these declarable at all: the bounding slot
        # carries the quantity, the bounded slot carries nothing
        slots = Dict(variable_slots(rel))
        @test slots[bd] === nothing
        @test slots[bg] === Q
        # ...so the relation is discoverable FROM the bounding quantity
        @test rel in AbstractQAtlas.relations_constraining(Q)
        @test quantities(rel) == (Q,)
    end
end

@testset "the bounded quantity has no type, so bounds_on finds nothing" begin
    # Not an oversight — the honest state of the vocabulary.  `bounds_on` matches
    # the BOUNDED slot's type; these have none, so the query is empty for all
    # eight while `relations_constraining` (checked above) is not.  When a
    # quantity for the measured side arrives, this flips, and it should.
    for (rel, _, _, _, Q) in _EIGHT
        @test isempty(bounds_on(Q))
    end
end

@testset "directions are the ones the literature states" begin
    # UPPER bounds: at the bound, slack 0; above it, violated; below it, satisfied.
    tsirelson = 2 * sqrt(2)
    @test slack(CHSHInequality(); S=tsirelson, S_max=tsirelson) == 0.0
    @test check(CHSHInequality(); S=2.0, S_max=tsirelson)             # classical S=2 obeys quantum
    @test !check(CHSHInequality(); S=3.0, S_max=tsirelson)            # above Tsirelson
    @test check(CHSHInequality(); S=3.0, S_max=4.0)                   # ...but fine no-signalling
    @test !check(MerminInequality(); M=4.5, M_max=4.0)
    @test check(MerminInequality(); M=4.0, M_max=4.0)                 # GHZ saturates
    @test !check(CloningFidelityBound(); F=0.9, F_max=5 // 6)
    @test check(CloningFidelityBound(); F=5 // 6, F_max=5 // 6)       # optimal cloner
    β = 2.0
    @test check(LyapunovChaosBound(); λ_L=1.0, λ_max=2π / β)
    @test !check(LyapunovChaosBound(); λ_L=2π / β + 0.1, λ_max=2π / β)
    @test check(BekensteinEntropyBound(); S=1.0, S_max=2π * 1.0 * 1.0)

    # LOWER bounds: the failure is being BELOW, which is the easy one to state
    # backwards — a protocol slower than the speed limit, or a key rate under
    # the proven achievable one.
    @test check(OrthogonalizationTimeBound(); τ=1.0, τ_min=0.5)
    @test !check(OrthogonalizationTimeBound(); τ=0.4, τ_min=0.5)
    @test slack(OrthogonalizationTimeBound(); τ=0.5, τ_min=0.5) == 0.0
    @test check(FastScramblingBound(); t_scr=2.0, t_min=1.0)
    @test !check(FastScramblingBound(); t_scr=0.5, t_min=1.0)
    @test check(SecretKeyRateBound(); r=0.6, r_min=0.5)
    @test !check(SecretKeyRateBound(); r=0.4, r_min=0.5)
end

@testset "solve returns the saturating value" begin
    @test solve(CHSHInequality(), Val(:S); S_max=2 * sqrt(2)) ≈ 2 * sqrt(2)
    @test solve(OrthogonalizationTimeBound(), Val(:τ); τ_min=0.5) ≈ 0.5
    # exactness survives: a Rational bound gives a Rational saturation
    @test solve(CloningFidelityBound(), Val(:F); F_max=5 // 6) === 5 // 6
end

@testset "an upper and a lower bound cannot be swapped without failing" begin
    # Guards the one thing a hand-written slack gets wrong silently.  Feed each
    # bound data that violates it, then feed the SAME data to a bound of the
    # opposite direction on the same slots — exactly one may pass.
    for (rel, bd, bg, dir, _) in _EIGHT
        violating = dir === :upper ? (2.0, 1.0) : (0.5, 1.0)
        satisfying = dir === :upper ? (0.5, 1.0) : (2.0, 1.0)
        @test !check(rel; NamedTuple{(bd, bg)}(violating)...)
        @test check(rel; NamedTuple{(bd, bg)}(satisfying)...)
    end
end
