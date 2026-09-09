using AbstractQAtlas
using ForwardDiff
using Test

@testset "peierls_current is −∂H/∂A by AD" begin
    A = VectorPotential(0.3)
    @test peierls_current(a -> a[1]^2, A) ≈ -0.6
    @test peierls_current(a -> cos(a[1]), VectorPotential(0.0)) ≈ 0.0
    @test peierls_current(a -> 3a[1], VectorPotential(1.0)) ≈ -3.0
    # The sign, isolated: a monotone increasing energy gives a NEGATIVE current, which is the
    # half of the convention that `|J|` can never see.
    @test peierls_current(a -> 3a[1], VectorPotential(1.0)) < 0
    @test peierls_current(a -> -3a[1], VectorPotential(1.0)) > 0
end

@testset "peierls_current agrees with a finite difference, and with the hand route" begin
    H(a) = 0.7 * cos(a[1] / 2) + 0.2 * sin(a[1]) - 0.1
    scalar(x) = H(VectorPotential(x))
    ε = 1.0e-6
    for a in (0.0, 0.35, -1.2)
        fd = (scalar(a + ε) - scalar(a - ε)) / (2ε)
        @test peierls_current(H, VectorPotential(a)) ≈ -fd atol = 1.0e-7
        # The two exported routes land on the same sign, so a consumer holding a derivative
        # already and one holding a function do not disagree.
        @test peierls_current(H, VectorPotential(a)) ≈
            current_from_hamiltonian_derivative(ForwardDiff.derivative(scalar, a))
    end
end

@testset "peierls_current carries no length unit of its own" begin
    # The seam differentiates with respect to whatever `A` it is handed; the unit is
    # `PeierlsConvention`'s business. Two models differing only by that unit give currents
    # differing by the same factor — the fact a consumer must not rediscover by fitting.
    Hsite(a) = cos(peierls_phase(PeierlsConvention(1), a))
    Hcell(a) = cos(peierls_phase(PeierlsConvention(2), a))
    A = 0.4
    @test peierls_current(Hcell, VectorPotential(A)) ≈
        peierls_current(Hsite, VectorPotential(A / 2)) / 2
    @test !isapprox(
        peierls_current(Hcell, VectorPotential(A)),
        peierls_current(Hsite, VectorPotential(A)),
    )
end

@testset "a higher dimension is still refused, with the backend loaded" begin
    # Both branches are exercised in this file — the extension's 1D method above and the
    # refusals here — so neither is left uncovered in whichever shard lacks ForwardDiff.
    msg = try
        peierls_current(a -> a[1] + a[2], VectorPotential(0.1, 0.2))
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("2-dimensional", msg)
    # And it is the dimension that refuses, not a missing backend: ForwardDiff IS loaded.
    @test !occursin("ForwardDiff", msg)
    @test peierls_current(a -> a[1]^2, VectorPotential(0.3)) ≈ -0.6
end
