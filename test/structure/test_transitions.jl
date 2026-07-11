# Phase-transition classification: the traits must encode the standard
# textbook distinctions (first vs second order vs BKT), each anchored in
# the Ehrenfest free-energy-derivative picture.

using AbstractQAtlas

@testset "transition type tree" begin
    for t in (FirstOrder(), ContinuousTransition(), KosterlitzThouless())
        @test t isa AbstractTransition
    end
end

@testset "Ehrenfest order (which free-energy derivative is singular)" begin
    @test ehrenfest_order(FirstOrder()) == 1
    @test ehrenfest_order(ContinuousTransition()) == 2
    @test ehrenfest_order(KosterlitzThouless()) == Inf
    @test ehrenfest_order(FirstOrder()) < ehrenfest_order(ContinuousTransition())
end

@testset "characterizing traits" begin
    # first order: latent heat, order parameter, NO critical exponents
    @test has_latent_heat(FirstOrder())
    @test has_order_parameter(FirstOrder())
    @test !has_critical_exponents(FirstOrder())

    # continuous: order parameter + critical exponents, no latent heat
    @test has_order_parameter(ContinuousTransition())
    @test has_critical_exponents(ContinuousTransition())
    @test !has_latent_heat(ContinuousTransition())

    # BKT: no local order parameter, no standard exponents, no latent heat
    @test !has_order_parameter(KosterlitzThouless())
    @test !has_critical_exponents(KosterlitzThouless())
    @test !has_latent_heat(KosterlitzThouless())

    # only the continuous transition carries the critical-exponent machinery
    withexp = filter(
        has_critical_exponents, [FirstOrder(), ContinuousTransition(), KosterlitzThouless()]
    )
    @test withexp == [ContinuousTransition()]
end
