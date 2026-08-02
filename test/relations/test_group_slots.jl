# Abstract GROUP slots under an explicit quantifier (design §6 / §8a).
#
# A parametric FAMILY needs no quantifier: `Susceptibility{(:x,:x)}` and
# `{(:z,:z)}` are one quantity at two indices, so a law written on the family is
# component-agnostic by construction and "check every component" is the only
# reading.  An abstract GROUP is different in kind — `MassGap`, `ChargeGap` and
# `SpinGap` are DIFFERENT quantities — and both readings are wanted:
#
#   EachOf{AbstractGap}   every gap is non-negative     — enumerate
#   AnyOf{AbstractGap}    THE gap sets ξ = 1/Δ          — the caller names one
#
# The measurement that forced the split is the `AnyOf` case below: applying the
# family rule to it reports the two irrelevant gaps as VIOLATED purely because
# they share a supertype.  So a BARE abstract slot is rejected and the author
# says which they mean.

using AbstractQAtlas
using Test
using AbstractQAtlas:
    AnyOf,
    EachOf,
    ambiguous_relations,
    applicable_relations,
    bag,
    check,
    quantities,
    relation_report,
    relations_constraining,
    residual,
    solve

# a downstream package's own group, which is the case §6 was actually about
module _GroupDownstream
using AbstractQAtlas
abstract type AbstractProbeQuantity <: AbstractQuantity end
struct ProbeAlpha <: AbstractProbeQuantity end
struct ProbeBeta <: AbstractProbeQuantity end
end
using ._GroupDownstream

function _unregister!(::Type{T}) where {T}
    return filter!(r -> !(r isa T), AbstractQAtlas._RELATION_REGISTRY)
end

@bound :test _EveryGapPositive(g::EachOf{AbstractGap} >= 0)
@relation :test _OneGapSetsXi(g::AnyOf{AbstractGap}, ξ) = g * ξ - 1
@relation :test _DownstreamGroupLaw(q::EachOf{_GroupDownstream.AbstractProbeQuantity}, x) =
    q - x

@testset "a BARE abstract group is refused, with both readings named" begin
    err = try
        @eval @relation :test _BareGroup(g::AbstractGap, x) = g - x
        nothing
    catch e
        sprint(showerror, e)
    end
    @test err !== nothing
    # the message must teach the fix, not just refuse
    @test occursin("EachOf", err)
    @test occursin("AnyOf", err)
    # a parametric family is still accepted with no quantifier — unchanged §8a
    @test_nowarn @eval @bound :test _FamilyStillFine(χ::Susceptibility >= 0)
    _unregister!(_FamilyStillFine)
end

@testset "EachOf enumerates: one row per member present" begin
    b = bag(MassGap => 0.5, ChargeGap => 1.5, SpinGap => -0.25)
    rows = [r for r in relation_report(b) if r.relation isa _EveryGapPositive]
    @test length(rows) == 3
    @test Set(r.subject.type for r in rows) == Set([MassGap, ChargeGap, SpinGap])
    # the negative one, and only it, fails
    @test [r.subject.type for r in rows if !r.pass] == [SpinGap]
    # a member absent from the bag produces no row (no fabricated subject)
    b2 = bag(MassGap => 0.5)
    rows2 = [r for r in relation_report(b2) if r.relation isa _EveryGapPositive]
    @test [r.subject.type for r in rows2] == [MassGap]
end

@testset "AnyOf is checkable but NOT auto-discovered" begin
    b = bag(MassGap => 0.5, ChargeGap => 2.0, SpinGap => 4.0)
    # the measurement that motivated the split: were this enumerated, ChargeGap and
    # SpinGap would be reported violated for sharing a supertype with MassGap
    @test !any(r -> r.relation isa _OneGapSetsXi, relation_report(b; ξ=2.0))
    @test !(_OneGapSetsXi() in applicable_relations(b; ξ=2.0))
    # ...and it is PENDING, not absent — the gap is listed
    @test _OneGapSetsXi() in ambiguous_relations(b; ξ=2.0)
    # naming the member makes it check, and it is the member that satisfies it
    @test check(_OneGapSetsXi(), b; subject=MassGap, ξ=2.0)
    @test !check(_OneGapSetsXi(), b; subject=ChargeGap, ξ=2.0)
    # a non-member is refused rather than silently matched
    @test_throws ErrorException residual(
        _OneGapSetsXi(), bag(Energy{:per_site} => 1.0); subject=Energy{:per_site}, ξ=2.0
    )
end

@testset "the two queries are disjoint, and neither is silently empty" begin
    b = bag(MassGap => 0.5, ChargeGap => 1.5)
    app = applicable_relations(b; ξ=2.0)
    amb = ambiguous_relations(b; ξ=2.0)
    @test !isempty(app)
    @test !isempty(amb)
    @test isempty(intersect(app, amb))
end

@testset "a quantified group enters the relation network" begin
    # the payoff: `quantities` unwraps the quantifier to the GROUP, so the reverse
    # index finds the relation from any concrete member
    @test quantities(_EveryGapPositive()) == (AbstractGap,)
    for M in (MassGap, ChargeGap, SpinGap)
        @test _EveryGapPositive() in relations_constraining(M)
    end
    # and a member of a DIFFERENT group does not match
    @test !(_EveryGapPositive() in relations_constraining(Energy{:per_site}))
end

@testset "a DOWNSTREAM group works the same — §6's actual ask" begin
    @test quantities(_DownstreamGroupLaw()) == (_GroupDownstream.AbstractProbeQuantity,)
    @test _DownstreamGroupLaw() in relations_constraining(_GroupDownstream.ProbeAlpha)
    b = bag(_GroupDownstream.ProbeAlpha => 3.0, _GroupDownstream.ProbeBeta => 5.0)
    rows = [r for r in relation_report(b; x=3.0) if r.relation isa _DownstreamGroupLaw]
    @test length(rows) == 2
    @test [r.pass for r in sort(rows; by=r -> string(r.subject.type))] == [true, false]
end

@testset "solve treats a group slot exactly as it treats a family slot" begin
    # §8a's decision, unchanged here: `solve`/`derive` targets stay concrete, and
    # a GENERIC slot is not solvable through by naming a member.  Measured on the
    # shipped family case first, so this is parity and not a new restriction —
    # `solve(SusceptibilityPositivity(), Susceptibility{(:z,:z)}, b)` throws the
    # same "has no variable of type" today.
    b = bag(MassGap => 0.5, ChargeGap => 1.5)
    bχ = bag(Susceptibility{(:z, :z)} => 0.7)
    @test_throws ErrorException solve(
        SusceptibilityPositivity(), Susceptibility{(:z, :z)}, bχ
    )
    @test_throws ErrorException solve(_EveryGapPositive(), MassGap, b)
    # a bare group as the TARGET is rejected with its own message, like a family
    @test_throws ErrorException solve(_EveryGapPositive(), AbstractGap, b)
    @test_throws ErrorException solve(SusceptibilityPositivity(), Susceptibility, bχ)
    # the symbol-keyed solve is unaffected — the slot is generic, the math is not
    @test solve(_EveryGapPositive(), Val(:g)) == 0
end

@testset "two generic slots are still refused (§8b not landed)" begin
    @test_throws Exception @eval @relation :test _TwoGeneric(
        a::EachOf{AbstractGap}, b::EachOf{AbstractVelocity}
    ) = a - b
end

foreach(_unregister!, (_EveryGapPositive, _OneGapSetsXi, _DownstreamGroupLaw))
