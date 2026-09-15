# The shipped examples, run and checked.
#
# An example nobody runs is worse than none: it reads as a promise and rots
# silently. This includes the script rather than restating it, so the two cannot
# disagree, and asserts against the EXACT values the physics fixes rather than
# against the digits the script happens to print.

using AbstractQAtlas
using Test: @test, @testset

include(joinpath(@__DIR__, "..", "..", "examples", "critical_entanglement.jl"))
using .CriticalEntanglementExample: block_entropy, central_charges

@testset "examples/critical_entanglement.jl" begin
    cs = central_charges()
    # The oracle is `c = 1` for a free-fermion chain, not a recorded output. The
    # gap is the finite-block correction at ℓ = 16, 64, which is why this is 1e-3.
    @test cs.nats ≈ 1 atol = 1e-3
    @test cs.bits_declared ≈ 1 atol = 1e-3
    # Undeclared bits are wrong by exactly the base, and the test says which
    # number it is rather than just "not 1": a wrong answer that happens to miss
    # by something else would otherwise pass for the same reason.
    @test cs.bits_undeclared ≈ 1 / log(2) atol = 1e-3
    @test !isapprox(cs.bits_undeclared, 1; atol=1e-3)
    # Declaring the convention has to change the answer, not merely be accepted.
    @test cs.bits_declared ≈ cs.nats atol = 1e-12

    # The entropies themselves are a measurement, so pin that they grow with the
    # block rather than pinning digits a BLAS version can move.
    @test block_entropy(64) > block_entropy(16) > block_entropy(8) > 0
end
