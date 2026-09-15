# Evaluating a response-genealogy edge along a chosen route.
#
# The routes here need no AD backend, so this file runs whether or not one is
# loaded; the AD/finite-difference agreement lives in test/ext/, which has one.

using AbstractQAtlas
using AbstractQAtlas: _route_order, _genealogy_derivative
using Test: @test, @test_throws, @testset

why(f) =
    try
        f()
        ""
    catch e
        sprint(showerror, e)
    end

F(h) = -log(2cosh(h))          # M = -F'(h) = tanh(h)
Φ(T) = -T * log(2cosh(1 / T))  # S = -Φ'(T)
kinked(x) = x < 0 ? x^2 : 2x^2 # f and f' continuous at 0, f'' jumps 2 -> 4

# A genealogy that roots somewhere other than a thermodynamic potential, so
# `_genealogy_derivative`'s root guard has something that can actually fire it.
struct RootProbeParent <: AbstractQuantity end
struct RootProbeQuantity <: AbstractQuantity end
function AbstractQAtlas.derivative_edge(::Type{RootProbeQuantity})
    return DerivativeEdge(RootProbeParent, Temperature)
end

# A route that breaks the contract: it reports a step and gives no way to change
# one. `observed_order` needs both, so this must surface rather than become a NaN.
struct BrokenStepRoute <: DerivativeRoute
    h::Float64
end
function AbstractQAtlas.nth_derivative(r::BrokenStepRoute, f, x, n::Integer)
    return nth_derivative(CentralDifference(r.h), f, x, n)
end
AbstractQAtlas.step_size(r::BrokenStepRoute) = r.h

# A route naming a backend that IS loaded, so the fallback must not blame the
# package for a method the route itself never defined.
struct LoadedBackendRoute <: DerivativeRoute end
AbstractQAtlas.backend_package(::LoadedBackendRoute) = :LinearAlgebra

@testset "a finite-difference route reaches the closed form" begin
    x = 0.3
    @test thermal_derivative(Magnetization(:z), F, x, CentralDifference(1e-2)) ≈ tanh(x) atol =
        1e-4
    @test thermal_derivative(Magnetization(:z), F, x, Richardson(1e-2)) ≈ tanh(x) atol =
        1e-11
    # Richardson is the reason to have it: same step, many orders closer.
    e_cd = abs(
        thermal_derivative(Magnetization(:z), F, x, CentralDifference(1e-2)) - tanh(x)
    )
    e_rd = abs(thermal_derivative(Magnetization(:z), F, x, Richardson(1e-2)) - tanh(x))
    @test e_rd < e_cd / 1e6
end

@testset "the sign and the order come from the genealogy, not from the route" begin
    # S = -∂F/∂T, one order in Temperature, and the minus sign of the tree.
    T = 1.7
    fd = (Φ(T + 1e-4) - Φ(T - 1e-4)) / 2e-4
    @test thermal_derivative(ThermalEntropy(), Φ, T, CentralDifference(1e-4)) ≈ -fd atol =
        1e-12
    # χ⁽²⁾ is a third derivative of F, and _route_order is what says so.
    @test _route_order(Susceptibility(:z, :z, :z)) == 3
    @test _route_order(Susceptibility(:z, :z)) == 2
    @test _route_order(Magnetization(:z)) == 1
    @test thermal_derivative(Susceptibility(:z, :z, :z), F, 0.3, Richardson(5e-2)) ≈
        -2 * sech(0.3)^2 * tanh(0.3) rtol = 1e-7
end

@testset "observed_order tells the truncation side from the roundoff side" begin
    # In the asymptotic regime a central difference shows its nominal order.
    @test observed_order(CentralDifference(1e-2), F, 0.3, 1) ≈ 2 atol = 0.05
    # A step small enough to be dominated by cancellation does not, even though
    # its error happens to be smaller. Written as "outside the band around 2"
    # rather than "below 2": within one decade of this step the quotient also
    # returns `Inf` (only the second difference vanishes) and `NaN` (both do), and
    # `!(o > 1.9)` is false for `Inf`, so that spelling would fail on a step 12%
    # away in log space with no platform difference needed.
    @test !(1.9 <= observed_order(CentralDifference(1e-9), F, 0.3, 1) <= 2.1)
    @test !(1.9 <= observed_order(CentralDifference(1e-11), F, 0.3, 1) <= 2.1)
    # The control the diagnostic needs: a function it should FAIL on. The second
    # derivative jumps at 0, which leaves the quotient first-order exactly: the
    # h^2 term of the central difference does not cancel, so D(h) = h/2.
    @test observed_order(CentralDifference(1e-3), kinked, 0.0, 1) ≈ 1 atol = 1e-9
    # AutoDiff reports no step, so there is nothing to halve. The message has to
    # say that, not the generic "no method".
    msg = try
        observed_order(AutoDiff(), F, 0.3, 1)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("step_size", msg)
    @test step_size(AutoDiff()) === nothing
    @test step_size(CentralDifference(1e-3)) == 1e-3
    @test step_size(Richardson(1e-2)) == 1e-2
    # Richardson carries a step too, so it must not fall out of `observed_order`
    # the way a closed Union over the routes that happened to exist would drop a
    # third one. Its value is noise once the route is at machine precision, so the
    # claim is that it RUNS and returns a number, not what the number is.
    @test isfinite(observed_order(Richardson(1e-1), F, 0.3, 1)) ||
        isnan(observed_order(Richardson(1e-1), F, 0.3, 1))
end

@testset "a route that carries a step must say so, not be named in a Union" begin
    # The contract is `step_size` + `with_step_size`, so a future route gets the
    # honest refusal instead of the false claim that it carries no step.
    @test with_step_size(CentralDifference(1e-2), 1e-3) == CentralDifference(1e-3)
    @test with_step_size(Richardson(1e-2; levels=4), 1e-3) == Richardson(1e-3; levels=4)
    @test_throws ErrorException with_step_size(AutoDiff(), 1e-3)
end

@testset "Richardson's levels is a knob, not a decoration" begin
    x, exact = 0.3, tanh(0.3)
    errs = [
        abs(
            thermal_derivative(Magnetization(:z), F, x, Richardson(1e-1; levels=L)) - exact
        ) for L in 2:5
    ]
    # STRICTLY falling: a `levels` that is ignored gives four equal errors, and
    # `issorted` counts ties as sorted, so it alone would pass that.
    @test all(errs[i] > errs[i + 1] for i in 1:(length(errs) - 1))
    @test errs[1] / errs[end] > 1e3
    # And a bad `levels` is refused rather than silently behaving as the default.
    @test_throws ArgumentError Richardson(1e-2; levels=1)
    @test_throws ArgumentError Richardson(1e-2; levels=0)
end

@testset "the root guard can fire" begin
    # Every shipped derivative_edge chains to FreeEnergy or GrandPotential, so
    # without a quantity rooted elsewhere this guard is unreachable and deleting it
    # changes nothing.
    @test potential_root(RootProbeQuantity()) === RootProbeParent
    msg = try
        thermal_derivative(RootProbeQuantity(), F, 0.3, CentralDifference(1e-3))
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("RootProbeParent", msg)
end

@testset "a route change cannot turn a refusal into a number" begin
    # The off-diagonal guard is on the AD path; it has to be on this one too.
    @test_throws ErrorException thermal_derivative(
        Susceptibility(:x, :y), F, 0.3, Richardson(1e-2)
    )
    # A quantity outside the genealogy has no edge to evaluate.
    @test_throws ErrorException thermal_derivative(
        FreeEnergy(), F, 0.3, CentralDifference(1e-3)
    )
    # Without a backend, AutoDiff refuses and names the routes that need none.
    msg = try
        thermal_derivative(Magnetization(:z), F, 0.3, AutoDiff())
        ""
    catch e
        sprint(showerror, e)
    end
    @test isempty(msg) || occursin("ForwardDiff", msg)
end

@testset "route arguments" begin
    @test_throws ArgumentError CentralDifference(0)
    @test_throws ArgumentError CentralDifference(-1e-3)
    @test_throws ArgumentError Richardson(1e-2; levels=1)
    @test_throws ArgumentError nth_derivative(CentralDifference(1e-3), F, 0.3, -1)
    for r in (CentralDifference(1e-3), Richardson(1e-2))
        @test nth_derivative(r, F, 0.3, 0) == F(0.3)
    end
    # `nth_derivative` is exported, so a caller can reach the AutoDiff route
    # directly rather than through `thermal_derivative`. Without a backend it has
    # to name the routes that need none, not fall through to the generic
    # "no method" of an unrecognised route.
    msg = try
        nth_derivative(AutoDiff(), F, 0.3, 2)
        ""
    catch e
        sprint(showerror, e)
    end
    @test isempty(msg) || occursin("CentralDifference", msg)
end

@testset "a report compares routes instead of trusting one" begin
    rows = derivative_report(
        Magnetization(:z), F, 0.3, (CentralDifference(1e-2), Richardson(1e-2), AutoDiff())
    )
    @test length(rows) == 3
    @test rows[1].order ≈ 2 atol = 0.05
    # A route that cannot run is a NaN row, not an aborted sweep: the usual cause
    # is a missing backend and the other rows are still the answer.
    ad = last(rows)
    @test isnan(ad.value) || isapprox(ad.value, tanh(0.3); atol=1e-12)
    @test isnan(ad.order)
    @test all(r -> isapprox(r.value, tanh(0.3); atol=1e-3), rows[1:2])
    # The row has to name the route it ran, or a report that always stored the
    # first route would read the same.
    @test [r.route for r in rows] == [CentralDifference(1e-2), Richardson(1e-2), AutoDiff()]
    # Only a missing backend is absorbed into a NaN row. A quantity with no
    # genealogy edge, and a guard the route itself raises, both propagate: a NaN
    # there would read as "install a package" for a mistake no package fixes.
    # Pinned by MESSAGE: `_route_order` has its own guard one line later, so a
    # type-only assertion passes whichever of the two fired.
    @test occursin(
        "is not a response function",
        why(
            () -> derivative_report(PartitionFunction(), F, 0.3, (CentralDifference(1e-2),))
        ),
    )
    @test_throws ErrorException derivative_report(
        Susceptibility(:x, :y), F, 0.3, (CentralDifference(1e-2),)
    )
end

@testset "the backend route's method belongs to the extension alone" begin
    # Defining `nth_derivative(::AutoDiff, ...)` in BOTH the package and the
    # extension is a method overwrite, which makes the extension fail to
    # precompile while every test here stays green, because the fallback path
    # still loads. So the package must own no method for that signature, and the
    # "which package" answer lives in a trait instead.
    owned = [
        m for m in methods(nth_derivative) if
        m.module === AbstractQAtlas && m.sig.parameters[2] === AutoDiff
    ]
    @test isempty(owned)
    @test backend_package(AutoDiff()) === :ForwardDiff
    @test backend_package(CentralDifference(1e-3)) === nothing
    @test backend_package(Richardson(1e-2)) === nothing
    # Without the extension the fallback has to say WHICH package, not "no method".
    msg = try
        nth_derivative(AutoDiff(), F, 0.3, 1)
        ""
    catch e
        sprint(showerror, e)
    end
    @test isempty(msg) || occursin("ForwardDiff", msg)
end

@testset "a route that breaks the step contract is not absorbed as a NaN" begin
    # The order column used to catch `ErrorException` as well as
    # `MissingRouteBackend`, so a route reporting a `step_size` with no
    # `with_step_size` produced the same NaN as one that legitimately has no step.
    # Different mistakes, and only one of them is the caller's.
    @test occursin(
        "with_step_size", why(() -> observed_order(BrokenStepRoute(1e-2), F, 0.3, 1))
    )
    @test_throws ErrorException derivative_report(
        Magnetization(:z), F, 0.3, (BrokenStepRoute(1e-2),)
    )
    # A route declaring no step is still a quiet NaN, the one case the column may
    # absorb, and its message says which half is missing.
    @test isnan(only(derivative_report(Magnetization(:z), F, 0.3, (AutoDiff(),))).order)
    @test occursin(
        "reports no `step_size`", why(() -> observed_order(AutoDiff(), F, 0.3, 1))
    )
end

@testset "a missing method is not blamed on a package that is loaded" begin
    # LinearAlgebra is loaded by this package, so reaching the fallback means the
    # ROUTE is incomplete. Telling its author to reinstall points away from that.
    msg = why(() -> nth_derivative(LoadedBackendRoute(), F, 0.3, 1))
    @test occursin("is loaded, but no method matched", msg)
    @test !occursin("which is not loaded", msg)
end

@testset "MissingRouteBackend cannot name an extension that does not exist" begin
    # It is exported, so an extension author can construct it. Built for a route
    # with no `backend_package` it used to render "needs the nothing extension".
    @test_throws ArgumentError MissingRouteBackend(CentralDifference(1e-3))
    @test_throws ArgumentError MissingRouteBackend(Richardson(1e-2))
    e = MissingRouteBackend(AutoDiff())
    @test e isa Exception
    @test occursin("MissingRouteBackend:", sprint(showerror, e))
    @test occursin("ForwardDiff", sprint(showerror, e))
end
