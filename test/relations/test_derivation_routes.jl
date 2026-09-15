# Every route to a target, rather than the first one the registry reaches.
#
# 66 of the 219 symbol-keyed outputs have more than one producing relation (38
# of 84 typed ones), up to nineteen, and `derive` runs whichever comes first.
#
# The trap this file exists to pin is not the disagreement, it is the vacuous
# agreement: derive the target first and the closure manufactures its own inputs
# from it, so several routes hand the supplied number back and read as
# confirmation of something nothing independent touched.

using AbstractQAtlas
using Test: @test, @test_throws, @testset

slope(c) = 2 * c / 6   # CFTEntanglementSlope: dS/dlnℓ = ncuts·c/6, ncuts = 2

@testset "every INDEPENDENT route appears, and only those" begin
    rows = derivation_routes(:c; dS_dlogℓ=slope(0.5), ncuts=2)
    # One measurement, one route. Letting the chain derive `c` first and then
    # rebuild the chord slope and the halved-chain difference FROM it returns
    # three rows that all read 0.5, which is this measurement counted three
    # times and would pass any agreement test put to it.
    @test length(rows) == 1
    @test nameof(typeof(first(rows).relation)) === :CFTEntanglementSlope
    @test first(rows).value ≈ 0.5
    # So asking for two independent routes here must fail: there is only one.
    @test_throws ErrorException derive_crosschecked(
        :c; dS_dlogℓ=slope(0.5), ncuts=2, min_routes=2
    )
end

@testset "the target is held out, so a route cannot confirm itself" begin
    # With only `c` and `ncuts`, the closure CAN manufacture the slopes: they are
    # derivable from the supplied `c`.
    reach = derivable(; c=0.9, ncuts=2)
    @test :dS_dlogℓ in reach
    @test :dS_dlogchord in reach
    # And yet no route is reported, because every one of them would be reading a
    # value built from the target. Three rows all returning 0.9 is what this
    # emptiness replaces.
    @test isempty(derivation_routes(:c; c=0.9, ncuts=2))
    # The same graph does reach `c` once something independent is supplied, so the
    # emptiness above is the holdout and not an unreachable target.
    @test !isempty(derivation_routes(:c; c=0.9, dS_dlogℓ=slope(0.5), ncuts=2))
end

@testset "a contradiction the first-route solver returns anyway" begin
    # `derive` hands back the supplied value without consulting anything.
    @test derive(:c; c=0.9, dS_dlogℓ=slope(0.5), ncuts=2) == 0.9
    msg = try
        derive_crosschecked(:c; c=0.9, dS_dlogℓ=slope(0.5), ncuts=2)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("disagree", msg)
    @test occursin("CFTEntanglementSlope", msg)   # the route is named
    @test occursin("supplied: 0.9", msg)          # and so is the value it contradicts
    # Consistent data passes, and returns the supplied value.
    @test derive_crosschecked(:c; c=0.5, dS_dlogℓ=slope(0.5), ncuts=2) == 0.5
    # No supplied target: the routes are the answer.
    @test derive_crosschecked(:c; dS_dlogℓ=slope(0.5), ncuts=2) ≈ 0.5
end

@testset "min_routes asks that a cross-check actually happened" begin
    # Default is permissive, and returns a number nothing checked.
    @test derive_crosschecked(:c; c=0.9, ncuts=2) == 0.9
    @test_throws ErrorException derive_crosschecked(:c; c=0.9, ncuts=2, min_routes=1)
    # It must not fire where a route does exist.
    @test derive_crosschecked(:c; c=0.5, dS_dlogℓ=slope(0.5), ncuts=2, min_routes=1) == 0.5
end

@testset "the typed door holds the target out the same way" begin
    β, Z, U = 0.8, 3.0, 0.4
    F = -log(Z) / β
    S = β * (U - F)                     # FreeEnergyLegendre: F = U - S/β
    # Two genuinely independent routes to F: through Z, and through the Legendre
    # transform. They agree here.
    good = bag(
        PartitionFunction => Z,
        InverseTemperature => β,
        Energy(:per_site) => U,
        ThermalEntropy => S,
    )
    rows = derivation_routes(FreeEnergy, good)
    @test length(rows) == 2
    @test Set(nameof(typeof(r.relation)) for r in rows) ==
        Set([:FreeEnergyFromZ, :FreeEnergyLegendre])
    @test all(r -> r.value ≈ F, rows)
    @test derive_crosschecked(FreeEnergy, good; min_routes=2) ≈ F

    # Break the entropy. `derive` still returns the RIGHT number, off the other
    # route, and says nothing about the input that contradicts it.
    bad = bag(
        PartitionFunction => Z,
        InverseTemperature => β,
        Energy(:per_site) => U,
        ThermalEntropy => S + 0.5,
    )
    @test derive(FreeEnergy, bad) ≈ F
    @test_throws ErrorException derive_crosschecked(FreeEnergy, bad)

    # The holdout, on the typed door. `FreeEnergyLegendre` produces BOTH `F` and
    # `S`, so with F supplied the closure can manufacture the S that gives F back.
    circular = bag(FreeEnergy => F, Energy(:per_site) => U, InverseTemperature => β)
    @test VariableKey(ThermalEntropy) in derivable(circular)   # it could be built
    @test isempty(derivation_routes(FreeEnergy, circular))     # and it is not used

    # β and T are one quantity under two names, so a bag carrying β must not let
    # `KelvinRelation` hand T back through a Peltier coefficient built from it.
    alias = bag(InverseTemperature => 0.5, Thermopower(:x, :x) => 3.0)
    @test VariableKey(PeltierCoefficient{(:x, :x)}) in derivable(alias)
    @test isempty(derivation_routes(Temperature, alias))
end
