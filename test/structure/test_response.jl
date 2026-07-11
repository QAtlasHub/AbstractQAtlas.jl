# The response-function derivative genealogy: the chain from any response
# back to the free energy, and the derivative order along it.
#
# Independent expectation: the textbook derivative structure —
# M = −∂F/∂h, χ = ∂M/∂h = ∂²F/∂h², U ⟵ F, C = ∂U/∂T — must be exactly
# what the genealogy reports, and every branch must root at FreeEnergy.

using AbstractQAtlas
using AbstractQAtlas:
    derivative_edge,
    differentiation_chain,
    potential_root,
    derivative_order,
    is_response,
    conjugate_field

@testset "genealogy edges" begin
    @test derivative_edge(MagnetizationZ()) == DerivativeEdge(FreeEnergy, MagneticField)
    @test derivative_edge(SusceptibilityZZ()) ==
        DerivativeEdge(MagnetizationZ, MagneticField)
    @test derivative_edge(SpecificHeat()) == DerivativeEdge(Energy, Temperature)
    @test derivative_edge(ThermalEntropy()) == DerivativeEdge(FreeEnergy, Temperature)
    # Energy is dispatched through its {G} parameter
    @test derivative_edge(Energy()) == DerivativeEdge(FreeEnergy, InverseTemperature)
    # roots / non-genealogy quantities have no edge
    @test derivative_edge(FreeEnergy()) === nothing
    @test derivative_edge(PartitionFunction()) === nothing
    @test derivative_edge(CorrelationLength()) === nothing
end

@testset "is_response" begin
    @test is_response(SusceptibilityZZ())
    @test is_response(MagnetizationX())
    @test is_response(SpecificHeat())
    @test !is_response(FreeEnergy())
    @test !is_response(PartitionFunction())
end

@testset "differentiation chain roots at the free energy" begin
    @test differentiation_chain(SusceptibilityZZ()) ==
        [SusceptibilityZZ, MagnetizationZ, FreeEnergy]
    @test differentiation_chain(MagnetizationX()) == [MagnetizationX, FreeEnergy]
    @test differentiation_chain(SpecificHeat()) == [SpecificHeat, Energy, FreeEnergy]
    @test differentiation_chain(FreeEnergy()) == [FreeEnergy]     # singleton root
    # every thermodynamic response roots at FreeEnergy
    for q in
        (SusceptibilityXX(), MagnetizationZ(), SpecificHeat(), ThermalEntropy(), Energy())
        @test potential_root(q) === FreeEnergy
    end
end

@testset "derivative order along the chain" begin
    # χ is a SECOND field-derivative of F; M is a first; C none
    @test derivative_order(SusceptibilityZZ(), MagneticField()) == 2
    @test derivative_order(MagnetizationZ(), MagneticField()) == 1
    @test derivative_order(SpecificHeat(), MagneticField()) == 0
    # C is a first T-derivative of U (one T-edge on the U→…→F chain)
    @test derivative_order(SpecificHeat(), Temperature()) == 1
    # F itself: zero of anything
    @test derivative_order(FreeEnergy(), MagneticField()) == 0
end

@testset "conjugate fields" begin
    @test conjugate_field(MagnetizationZ()) == MagneticField()
    @test conjugate_field(MagnetizationX()) == MagneticField()
    @test conjugate_field(ThermalEntropy()) == Temperature()
end

@testset "exact-formula companions (the relations the edges point at)" begin
    # χ = ∂M/∂h  (SusceptibilityResponse) and χ = β·Var(M) (SusceptibilityFDT)
    # are the same response two ways — both live in the registry now.
    @test check(SusceptibilityResponse(); χ=2.5, dM_dh=2.5)
    @test check(MagnetizationResponse(); M=1.0, dF_dh=-1.0)   # M = −∂F/∂h
    # the definitional and statistical routes agree on a shared value:
    dM_dh = 0.8               # ∂M/∂h at some point
    var_M = 1.6
    β = 0.5                   # β·Var(M) = 0.8 = ∂M/∂h  ⇒ same χ
    χ_def = solve(SusceptibilityResponse(), Val(:χ); dM_dh=dM_dh)
    χ_stat = solve(SusceptibilityFDT(), Val(:χ); var_M=var_M, β=β, N=1)
    @test χ_def ≈ χ_stat
end
