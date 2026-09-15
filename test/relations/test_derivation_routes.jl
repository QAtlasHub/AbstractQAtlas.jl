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
using AbstractQAtlas: _disagreement, derivation_steps, typed_derivation_steps

why(f) =
    try
        f()
        ""
    catch e
        sprint(showerror, e)
    end
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
    @test occursin("differ by 0.4", msg)      # the size, not just that it threw
    @test occursin("CFTEntanglementSlope", msg)   # the route is named
    @test occursin("supplied: 0.9", msg)          # and so is the value it contradicts
    # Consistent data passes, and returns the supplied value.
    @test derive_crosschecked(:c; c=0.5, dS_dlogℓ=slope(0.5), ncuts=2) == 0.5
    # No supplied target: the routes are the answer.
    @test derive_crosschecked(:c; dS_dlogℓ=slope(0.5), ncuts=2) ≈ 0.5
end

@testset "min_routes defaults to demanding that a cross-check happened" begin
    # Data affording no independent route is refused by default: returning 0.9 from
    # it is what the verb's name would otherwise be claiming it had checked.
    @test_throws ErrorException derive_crosschecked(:c; c=0.9, ncuts=2)
    @test derive_crosschecked(:c; c=0.9, ncuts=2, min_routes=0) == 0.9   # opt-out
    # One route is enough to compare a supplied value against.
    @test derive_crosschecked(:c; c=0.5, dS_dlogℓ=slope(0.5), ncuts=2) == 0.5
    # And two can be demanded where one is not evidence.
    @test_throws ErrorException derive_crosschecked(
        :c; c=0.5, dS_dlogℓ=slope(0.5), ncuts=2, min_routes=2
    )
end

@testset "a route that raised is reported, not dropped" begin
    # Z must be positive: it is a sum of Boltzmann weights. `FreeEnergyFromZ` needs
    # log(Z) and throws, which is the relation being PREVENTED from disagreeing.
    # Dropping it silently leaves one route and a clean "cross-checked" answer.
    impossible = bag(
        PartitionFunction => -2.0,
        InverseTemperature => 1.0,
        Energy(:per_site) => -0.4,
        ThermalEntropy => 0.3,
    )
    rows = derivation_routes(FreeEnergy, impossible)
    @test length(rows) == 2
    threw = only(r for r in rows if r.error !== nothing)
    @test nameof(typeof(threw.relation)) === :FreeEnergyFromZ
    @test occursin("DomainError", threw.error)
    @test threw.value === nothing
    # The row DISPLAYS as a failure. Printing `r.value` unconditionally would show
    # `nothing`, and the refusal message is built from `string.(rows)`.
    @test occursin("THREW", string(threw))
    @test occursin("DomainError", string(threw))
    # Neither state and both states are unrepresentable, so no consumer's
    # `error === nothing` branch can be wrong about what `value` holds.
    @test_throws ArgumentError DerivationRouteRow(threw.relation, Any[], nothing, nothing)
    @test_throws ArgumentError DerivationRouteRow(threw.relation, Any[], 1.0, "boom")
    msg = try
        derive_crosschecked(FreeEnergy, impossible)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("never got to disagree", msg)
    @test occursin("FreeEnergyFromZ", msg)
    # A relation merely declining to be solved for a slot is NOT a broken route, or
    # every ordinary call would refuse. The consistent bag still passes.
    β, Z, U = 0.8, 3.0, 0.4
    F = -log(Z) / β
    ok = bag(
        PartitionFunction => Z,
        InverseTemperature => β,
        Energy(:per_site) => U,
        ThermalEntropy => β * (U - F),
    )
    @test all(r -> r.error === nothing, derivation_routes(FreeEnergy, ok))
    @test derive_crosschecked(FreeEnergy, ok; min_routes=2) ≈ F
end

@testset "a NaN route is a degeneration, not an agreement" begin
    # Every difference against NaN is NaN, so an `isnan` short-circuit meant for
    # "fewer than two values" would read a degenerate route as agreement.
    msg = try
        derive_crosschecked(:c; c=NaN, dS_dlogℓ=slope(0.5), ncuts=2)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("NaN", msg)
    @test occursin("nothing was compared", msg)
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

@testset "the tolerance is isapprox's, not a floor that goes absolute below one" begin
    # Two routes returning +1e-12 and -1e-12 is a sign flip. Dividing the
    # difference by `max(maximum(abs, vs), 1)` reports it as 2e-12 and passes it
    # at any sane rtol, which is why the comparison is `d <= atol + rtol*m`.
    d, m = _disagreement([1e-12, -1e-12])
    @test (d, m) == (2e-12, 1e-12)
    @test !(d <= 0 + 1e-8 * m)
    # A genuine agreement to 1e-9 relative still passes.
    d2, m2 = _disagreement([1.0, 1.0 + 1e-9])
    @test d2 <= 0 + 1e-8 * m2
    # Fewer than two values is not agreement. It gets `nothing`, not a NaN that a
    # degenerate route would also produce.
    @test _disagreement([1.0]) === nothing
    @test _disagreement(Any[]) === nothing

    # `atol` is how a caller says their routes are noise-dominated, and it is the
    # only thing that lets the broken-entropy bag through.
    β, Z, U = 0.8, 3.0, 0.4
    F = -log(Z) / β
    bad = bag(
        PartitionFunction => Z,
        InverseTemperature => β,
        Energy(:per_site) => U,
        ThermalEntropy => β * (U - F) + 0.5,
    )
    @test_throws ErrorException derive_crosschecked(FreeEnergy, bad)
    @test derive_crosschecked(FreeEnergy, bad; atol=1.0) ≈ F
end

@testset "multiple producers is the common case the section comment claims" begin
    # The comment above `derivation_routes` carries measured counts. Pinned as
    # floors rather than exact numbers: adding a relation that produces an
    # already-produced output raises them, and a floor does not rot for that.
    function multi(steps)
        d = Dict{Any,Set{Any}}()
        for st in steps
            push!(get!(d, st.output, Set{Any}()), nameof(typeof(st.relation)))
        end
        return count(v -> length(v) > 1, values(d)), maximum(length, values(d))
    end
    n_sym, max_sym = multi(derivation_steps())
    n_typed, max_typed = multi(typed_derivation_steps())
    @test n_sym >= 60
    @test max_sym >= 15
    @test n_typed >= 35
    @test max_typed >= 12
end

@testset "the declining/broken split holds across the whole registry" begin
    # `_route_declined` separates "this relation cannot be applied here" (skip)
    # from "it applied and the data broke it" (report). Getting it wrong in either
    # direction is invisible: too narrow and a broken input is silently dropped,
    # too wide and ordinary calls refuse.
    #
    # Both framework declination shapes, on fixtures that actually produce them.
    @test isempty(derivation_routes(:β; C=1.0, var_E=1.0))          # solve: not affine
    @test isempty(                                                   # untyped slot absent
        derivation_routes(
            Temperature, bag(InverseTemperature => 0.5, Thermopower(:x, :x) => 3.0)
        ),
    )
    # A relation's own physics guard is about the DATA and must be reported, not
    # skipped. `ncuts = 0` makes the residual independent of `c`.
    ncz = derivation_routes(:c; dS_dlogℓ=0.1667, ncuts=0)
    @test length(ncz) == 1
    @test occursin("ncuts = 0", only(ncz).error)

    # The sweep: hand every target's inputs the same nonsense value and check that
    # nothing the framework MEANT as a declination ends up in the reported set. A
    # seventh declination site added to interface.jl and not registered here fails
    # this, which is the whole point of the assertion.
    leaked = String[]
    total = 0
    for t in unique(s.output for s in derivation_steps())
        ins = unique(
            vcat([collect(s.inputs) for s in derivation_steps() if s.output === t]...)
        )
        rows = try
            derivation_routes(t; (i => 0.7 for i in ins)...)
        catch
            continue
        end
        total += length(rows)
        for r in rows
            r.error === nothing && continue
            (startswith(r.error, "solve:") || occursin("(untyped slot)", r.error)) &&
                push!(leaked, r.error)
        end
    end
    @test total > 300               # the sweep reached the registry, not two rows
    @test isempty(leaked)
end

@testset "a guard on the solver's probe must not make a slot unreachable" begin
    # `solve` probes its target at 0, 1, 2. A relation guarding one of those points
    # then refuses for EVERY input, and the refusal names a value the caller never
    # supplied. Both shapes below did that, and the fix is the closed form each
    # relation already had in prose.
    #
    # `EntanglementSpectrumCorrelation` is `ε - log((1-ζ)/ζ)`: at the probe ζ = 2
    # the argument of `log` is -0.5. Its docstring carried `ζ = 1/(e^ε + 1)` already.
    @test derive(:ζ; ε=0.7) ≈ 1 / (exp(0.7) + 1)
    @test derive_crosschecked(:ζ; ε=0.7) ≈ 1 / (exp(0.7) + 1)
    @test isempty([r for r in derivation_routes(:ζ; ε=0.7) if r.error !== nothing])

    # The slope relations guard `ncuts = 0`, which is exactly the first probe.
    @test derive(:ncuts; dS_dlogℓ=slope(0.5), c=0.5) ≈ 2
    @test derive_crosschecked(:ncuts; dS_dlogℓ=slope(0.5), c=0.5) ≈ 2
    @test isempty([
        r for
        r in derivation_routes(:ncuts; dS_dlogℓ=slope(0.5), c=0.5) if r.error !== nothing
    ])
    # The guard still holds, now read off the ANSWER rather than the probe, and it
    # reaches the caller. Through `derive` it does not: that verb drops the route
    # and reports the target unreachable, which is the difference this layer makes.
    @test occursin("ncuts = 0", only(derivation_routes(:ncuts; dS_dlogℓ=0.0, c=0.5)).error)
    @test occursin("ncuts = 0", why(() -> derive_crosschecked(:ncuts; dS_dlogℓ=0.0, c=0.5)))
    @test occursin("c = 0", why(() -> derive_crosschecked(:ncuts; dS_dlogℓ=0.3, c=0.0)))
    @test occursin("not reachable", why(() -> derive(:ncuts; dS_dlogℓ=0.0, c=0.5)))
end
