# Evaluating a response-genealogy edge along a chosen route.
#
# The routes here need no AD backend, so this file runs whether or not one is
# loaded; the AD/finite-difference agreement lives in test/ext/, which has one.

using AbstractQAtlas
using AbstractQAtlas: _route_order, _genealogy_derivative
using Test: @test, @test_throws, @testset

F(h) = -log(2cosh(h))          # M = -F'(h) = tanh(h)
Φ(T) = -T * log(2cosh(1 / T))  # S = -Φ'(T)
kinked(x) = x < 0 ? x^2 : 2x^2 # derivative discontinuous at 0

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
    # its error happens to be smaller. Stated as "does not show the nominal
    # order" rather than a number, because the value there is roundoff and would
    # be a different number on another machine (NaN included, when the successive
    # differences both vanish).
    @test !(observed_order(CentralDifference(1e-9), F, 0.3, 1) > 1.9)
    # The control the diagnostic needs: a function it should FAIL on. A kink at
    # the evaluation point leaves the quotient first-order, exactly.
    @test observed_order(CentralDifference(1e-3), kinked, 0.0, 1) ≈ 1 atol = 1e-9
    @test_throws ErrorException observed_order(AutoDiff(), F, 0.3, 1)
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
end
