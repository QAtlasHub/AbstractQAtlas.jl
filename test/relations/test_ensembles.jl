# Ensemble relations vs INDEPENDENT exact / statistical constructions.
#
# MicrocanonicalTemperature is checked through ensemble equivalence: for a
# solvable model the microcanonical ∂S/∂E must reproduce the canonical β
# at the matching energy.  CanonicalTPQ is checked both as an exact
# identity and by the random-state average it estimates.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve
using LinearAlgebra, Random

@testset "MicrocanonicalTemperature: β = ∂S/∂E and ensemble equivalence" begin
    # N two-level systems (energies 0 / ε): binary entropy per site
    # s(x) = −x ln x − (1−x) ln(1−x), energy per site e = ε x.  The
    # microcanonical β = ∂s/∂e must equal the canonical β at the canonical
    # occupation x = 1/(e^{βε}+1) — that identity IS ensemble equivalence.
    ε = 1.0
    s(x) = (x <= 0 || x >= 1) ? 0.0 : -x * log(x) - (1 - x) * log(1 - x)
    for β in (0.3, 1.0, 2.0, 3.5)
        x = 1 / (exp(β * ε) + 1)                     # canonical occupation
        h = 1e-6
        dS_dE = (s(x + h) - s(x - h)) / (2 * (ε * h))  # ∂s/∂e, finite difference
        @test check(MicrocanonicalTemperature(); β=β, dS_dE=dS_dE, atol=1e-4)
        # and solve recovers β from the microcanonical derivative
        @test solve(MicrocanonicalTemperature(), Val(:β); dS_dE=dS_dE) ≈ β atol = 1e-4
    end
end

@testset "CanonicalTPQ: Z = D·⟨ψ₀|e^{−βH}|ψ₀⟩" begin
    # exact form: the ideal random-state weight is Z/D, so the relation closes
    levels = [0.0, 0.7, 1.3, 2.1, 2.9]
    D = length(levels)
    β = 0.9
    Z = sum(exp(-β * E) for E in levels)
    @test check(CanonicalTPQ(); Z=Z, tpq_weight=Z / D, D=D, atol=1e-13)
    @test solve(CanonicalTPQ(), Val(:Z); tpq_weight=Z / D, D=D) ≈ Z
    @test solve(CanonicalTPQ(), Val(:tpq_weight); Z=Z, D=D) ≈ Z / D

    # statistical: the random-state average of ⟨ψ₀|e^{−βH}|ψ₀⟩ converges to
    # Z/D, so D·average estimates Z.  Diagonal H ⇒ weight = Σ|c_i|² e^{−βE_i}.
    rng = MersenneTwister(12345)
    w = exp.(-β .* levels)
    M = 20000
    acc = 0.0
    for _ in 1:M
        c = randn(rng, ComplexF64, D)
        c ./= norm(c)                       # Haar-random normalized state
        acc += real(sum(abs2.(c) .* w))
    end
    Z_est = D * acc / M
    @test isapprox(Z_est, Z; rtol=0.03)     # 20k samples ⇒ a few % is comfortable
end

@testset "type-keyed: ensemble" begin
    @test quantities(CanonicalTPQ()) == (PartitionFunction,)
    # MicrocanonicalTemperature keeps its S↔E graph edges via also_constrains
    @test Set(quantities(MicrocanonicalTemperature())) == Set((ThermalEntropy, Energy))
    # CanonicalTPQ: Z = D·⟨tpq_weight⟩ via bag (Z by type, weight + dimension via extras)
    @test check(
        CanonicalTPQ(), bag(PartitionFunction => 6.0); tpq_weight=2.0, D=3.0, atol=1e-12
    )
    # MicrocanonicalTemperature β = ∂S/∂E via bag (β field; dS_dE supplied)
    @test check(
        MicrocanonicalTemperature(), bag(InverseTemperature => 0.5); dS_dE=0.5, atol=1e-12
    )
end
