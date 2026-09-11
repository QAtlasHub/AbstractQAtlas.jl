# Scaling relations vs INDEPENDENT exact exponent sets.
#
# The 2D Ising rationals and the mean-field set are exact, independently
# known values (Onsager/Yang lattice solutions; Landau theory) — the
# relations must hold with residual ≡ 0 in EXACT arithmetic, which also
# pins the no-float-promotion contract.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve
using LinearAlgebra

const ISING2D = (α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4)
const MEANFIELD = (α=0//1, β=1//2, γ=1//1, δ=3//1, ν=1//2, η=0//1)

@testset "2D Ising exact rationals: residuals ≡ 0//1, types preserved" begin
    r1 = residual(Rushbrooke(); α=ISING2D.α, β=ISING2D.β, γ=ISING2D.γ)
    r2 = residual(Widom(); β=ISING2D.β, γ=ISING2D.γ, δ=ISING2D.δ)
    r3 = residual(Fisher(); γ=ISING2D.γ, ν=ISING2D.ν, η=ISING2D.η)
    r4 = residual(Josephson(); α=ISING2D.α, ν=ISING2D.ν, d=2)
    for r in (r1, r2, r3, r4)
        @test r isa Rational       # exact-arithmetic contract
        @test r == 0//1            # exact, not ≈
    end
    @test exponents_consistent(ISING2D; d=2)          # atol = 0: exact gate
    @test all(iszero, values(exponent_residuals(ISING2D; d=2)))
end

@testset "mean-field set at the upper critical dimension d = 4" begin
    @test exponents_consistent(MEANFIELD; d=4)
    # hyperscaling FAILS off the upper critical dimension — that failure
    # is physics (mean-field violates Josephson for d < 4), so the gate
    # must flag it:
    @test !check(Josephson(); α=MEANFIELD.α, ν=MEANFIELD.ν, d=3)
end

@testset "3D Ising bootstrap values within quoted precision" begin
    # Kos–Poland–Simmons-Duffin–Vichi (2016) determine (Δ_σ, Δ_ε); the
    # standard exponent set derives from them, so this is a consistency
    # check of our arithmetic against the published rounded values (the
    # genuinely independent exactness test is the 2D rational case above).
    nt3 = (α=0.11009, β=0.32642, γ=1.23708, δ=4.78984, ν=0.62999, η=0.03631)
    # residuals limited by the 5-digit rounding of the published values:
    @test exponents_consistent(nt3; d=3, atol=2e-4)
end

@testset "solve: every solvable variable round-trips exactly" begin
    @test solve(Rushbrooke(), Val(:α); β=1//8, γ=7//4) == 0//1
    @test solve(Rushbrooke(), Val(:β); α=0//1, γ=7//4) == 1//8
    @test solve(Rushbrooke(), Val(:γ); α=0//1, β=1//8) == 7//4
    @test solve(Widom(), Val(:β); γ=7//4, δ=15//1) == 1//8
    @test solve(Widom(), Val(:γ); β=1//8, δ=15//1) == 7//4
    @test solve(Widom(), Val(:δ); β=1//8, γ=7//4) == 15//1
    @test solve(Fisher(), Val(:γ); ν=1//1, η=1//4) == 7//4
    @test solve(Fisher(), Val(:ν); γ=7//4, η=1//4) == 1//1
    @test solve(Fisher(), Val(:η); γ=7//4, ν=1//1) == 1//4
    @test solve(Josephson(), Val(:α); ν=1//1, d=2) == 0//1
    @test solve(Josephson(), Val(:ν); α=0//1, d=2) == 1//1
    @test solve(Josephson(), Val(:d); α=0//1, ν=1//1) == 2//1
    # solve-then-residual ≡ 0 (definitional round trip)
    γ = solve(Widom(), Val(:γ); β=1//8, δ=15//1)
    @test residual(Widom(); β=1//8, γ=γ, δ=15//1) == 0//1
end

@testset "a wrong exponent set fails the gate" begin
    bad = (α=0//1, β=1//8, γ=3//2, δ=15//1, ν=1//1, η=1//4)   # γ ≠ 7/4
    @test !exponents_consistent(bad; d=2)
    @test !check(Widom(); β=bad.β, γ=bad.γ, δ=bad.δ)
end

@testset "activated scaling is a different law, not another value of z" begin
    # The gap of the CLEAN critical Ising chain, exactly: the Majorana matrix
    # A of H = -Σ Z_j Z_{j+1} - Σ X_j has smallest singular value E₁ - E₀.
    # Deterministic — no sampling anywhere in this testset.
    function gap(N)
        A = zeros(2N, 2N)
        for j in 1:N
            A[2j - 1, 2j] = -2.0
        end
        for j in 1:(N - 1)
            A[2j, 2j + 1] = -2.0
        end
        return minimum(svdvals(A - transpose(A)))
    end
    localslope(y, x, k) = (y[k + 1] - y[k]) / (x[k + 1] - x[k])

    Ns = (16, 32, 64, 128, 256)
    lnN = log.(collect(Float64, Ns))
    Δ = [gap(N) for N in Ns]
    z_eff = [-localslope(log.(Δ), lnN, k) for k in 1:4]
    ψ_eff = [localslope(log.(-log.(Δ)), lnN, k) for k in 1:4]

    # At a conventional critical point a constant `z` EXISTS and a constant `ψ`
    # does not.  Measured, deterministic:
    #
    #   z_eff  0.9776  0.9888  0.9944  0.9972   → spread 1.020, tending to 1
    #   ψ_eff  0.4941  0.3711  0.2964  0.2464   → spread 2.005, still falling
    #
    # This is the whole reason the two relations cannot be interchanged, and it
    # needs no disorder average to see.
    @test maximum(z_eff) / minimum(z_eff) < 1.05
    @test maximum(ψ_eff) / minimum(ψ_eff) > 1.8
    @test check(DynamicalScaling(); dlogΔ_dlogξ=(-z_eff[end]), z=1.0, atol=0.005)

    # The converse direction, on a family that IS activated by construction:
    # Δ = exp(-a ξ^ψ).  Here the activated slope is flat and `z_eff` is the one
    # that runs away, so neither law can stand in for the other.
    a, ψtrue = 0.7, 0.5
    ξ = collect(Float64, Ns)
    Δa = exp.(-a .* ξ .^ ψtrue)
    ψa = [localslope(log.(-log.(Δa)), lnN, k) for k in 1:4]
    za = [-localslope(log.(Δa), lnN, k) for k in 1:4]
    @test all(x -> isapprox(x, ψtrue; atol=1e-12), ψa)      # exact: ψ is the slope
    # No constant z exists, and the way it runs away is derived rather than
    # eyeballed: `−d(ln Δ)/d(ln ξ) = ψ ln(1/Δ) = ψ a ξ^ψ`, so each DOUBLING of
    # the size multiplies the effective z by exactly `2^ψ`.
    @test all(k -> za[k + 1] / za[k] ≈ 2^ψtrue, 1:3)
    @test check(ActivatedDynamicalScaling(); dloglogΔ_dlogξ=ψa[end], ψ=ψtrue, atol=1e-12)
    @test solve(ActivatedDynamicalScaling(), Val(:ψ); dloglogΔ_dlogξ=ψa[end]) ≈ ψtrue

    # And each law REJECTS the other's data at the far end of the range.
    @test !check(DynamicalScaling(); dlogΔ_dlogξ=(-za[end]), z=za[1], atol=0.5)
    @test !check(
        ActivatedDynamicalScaling(); dloglogΔ_dlogξ=ψ_eff[end], ψ=ψ_eff[1], atol=0.05
    )

    # Both exponents are typed subjects, so the registry can find the relations.
    @test variable_types(ActivatedDynamicalScaling()) == (ActivatedExponent,)
    @test variable_types(DynamicalScaling()) == (DynamicalExponent,)
    @test ActivatedDynamicalScaling() in relations_constraining(ActivatedExponent)
    @test DynamicalScaling() in relations_constraining(DynamicalExponent)
end

@testset "both dynamical laws are reachable from the gap they are about" begin
    # `Δ` enters each through the SUPPLIED derivative and so has no slot; without
    # the `also_constrains` link, asking what constrains a gap would not mention
    # how it closes — and would not offer the activated alternative at all.
    for rel in (DynamicalScaling(), ActivatedDynamicalScaling())
        @test MassGap in quantities(rel)
        @test rel in relations_constraining(MassGap)
    end
    @test DynamicalExponent in quantities(DynamicalScaling())
    @test ActivatedExponent in quantities(ActivatedDynamicalScaling())
end
