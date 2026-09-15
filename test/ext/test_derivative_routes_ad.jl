# The two routes on the same genealogy edge. `_genealogy_derivative` is shared,
# so this is what would catch the sign or the order drifting between them.

using AbstractQAtlas
using ForwardDiff: ForwardDiff
using Test: @test, @testset

Fad(h) = -log(2cosh(h))
Φad(T) = -T * log(2cosh(1 / T))
Uad(T) = 1.5T^2
βFad(β) = -log(2cosh(β)) / β
Ωad(μ) = -log(1 + exp(μ))

@testset "AutoDiff as a route is the backend method" begin
    @test thermal_derivative(Magnetization(:z), Fad, 0.3, AutoDiff()) ==
        thermal_derivative(Magnetization(:z), Fad, 0.3)
    @test thermal_derivative(SpecificHeat(), Uad, 2.0, AutoDiff()) ==
        thermal_derivative(SpecificHeat(), Uad, 2.0)
    @test thermal_derivative(Energy(), βFad, 0.7, AutoDiff()) ==
        thermal_derivative(Energy(), βFad, 0.7)
end

@testset "Richardson agrees with AD across the genealogy" begin
    cases = (
        (Magnetization(:z), Fad, 0.3, Richardson(1e-2)),
        (Susceptibility(:z, :z), Fad, 0.3, Richardson(2e-2)),
        (Susceptibility(:z, :z, :z), Fad, 0.3, Richardson(5e-2)),
        (ThermalEntropy(), Φad, 1.7, Richardson(1e-2)),
        (ParticleNumber(), Ωad, 0.4, Richardson(1e-2)),
        (SpecificHeat(), Uad, 2.0, Richardson(1e-2)),
        (Energy(), βFad, 0.7, Richardson(1e-2)),
    )
    for (q, pot, x, route) in cases
        ad = thermal_derivative(q, pot, x, AutoDiff())
        fd = thermal_derivative(q, pot, x, route)
        # Sign first: a shared genealogy is the claim, and a sign flip would still
        # pass a loose magnitude comparison on a symmetric point.
        @test sign(ad) == sign(fd)
        @test isapprox(fd, ad; rtol=1e-6)
    end
end

@testset "the report carries the backend row exactly" begin
    rows = derivative_report(
        Magnetization(:z), Fad, 0.3, (AutoDiff(), Richardson(1e-2), CentralDifference(1e-2))
    )
    @test rows[1].value == thermal_derivative(Magnetization(:z), Fad, 0.3)
    @test isnan(rows[1].order)          # no step to halve
    @test rows[3].order ≈ 2 atol = 0.05
    @test maximum(r -> abs(r.value - rows[1].value), rows) < 1e-4
end
