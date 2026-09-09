using AbstractQAtlas
using ForwardDiff
using Test

@testset "peierls_current is −∂H/∂A by AD" begin
    @test peierls_current(a -> a^2, 0.3) ≈ -0.6
    @test peierls_current(cos, 0.0) ≈ 0.0
    @test peierls_current(a -> 3a, 1.0) ≈ -3.0
    # The sign, isolated: a monotone increasing energy gives a NEGATIVE current, which is the
    # half of the convention that `|J|` can never see.
    @test peierls_current(a -> 3a, 1.0) < 0
    @test peierls_current(a -> -3a, 1.0) > 0
end

@testset "peierls_current agrees with a finite difference, and with the hand route" begin
    # Differentiating at the seam has to give what a caller would get by differencing their
    # own energy — that is the claim which makes hand-coding unnecessary.
    H(a) = 0.7 * cos(a / 2) + 0.2 * sin(a) - 0.1
    ε = 1.0e-6
    for A in (0.0, 0.35, -1.2)
        fd = (H(A + ε) - H(A - ε)) / (2ε)
        @test peierls_current(H, A) ≈ -fd atol = 1.0e-7
        # …and the two exported routes agree, so a consumer holding a derivative already and
        # one holding a function land on the same sign.
        @test peierls_current(H, A) ≈
            current_from_hamiltonian_derivative(ForwardDiff.derivative(H, A))
    end
end

@testset "peierls_current carries no length unit of its own" begin
    # The seam differentiates with respect to whatever `A` the caller passes; the unit is
    # `PeierlsConvention`'s business. Two models differing only by that unit therefore give
    # currents differing by the same factor — which is the fact a consumer must not
    # rediscover by fitting.
    Hsite(a) = cos(peierls_phase(PeierlsConvention(1), a))
    Hcell(a) = cos(peierls_phase(PeierlsConvention(2), a))
    A = 0.4
    @test peierls_current(Hcell, A) ≈ peierls_current(Hsite, A / 2) / 2
    # Not equal, which is the whole point of naming the unit.
    @test !isapprox(peierls_current(Hcell, A), peierls_current(Hsite, A))
end

@testset "the fallback still refuses, with the extension loaded" begin
    # Both branches of `peierls_current` are exercised in this file: the extension's method
    # above, and the fallback here through `invoke`. Reaching only one leaves the other
    # uncovered in whichever shard happens not to load ForwardDiff.
    @test_throws ErrorException invoke(peierls_current, Tuple{Any,Any}, a -> a^2, 0.3)
    # And the extension's method is genuinely a different one, not the fallback rethrown.
    @test peierls_current(a -> a^2, 0.3) ≈ -0.6
end
