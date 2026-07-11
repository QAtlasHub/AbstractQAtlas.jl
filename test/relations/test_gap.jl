# Gap ↔ correlation length, and the dynamical exponent.
#
# ξ = v/Δ is verified against the real-space decay of a relativistic
# dispersion; the dynamical exponent is read off a synthetic Δ(ξ) power
# law.  This file also pins the generic-solve guard against a probe that
# blows up (ξ = v/Δ is non-affine in Δ).

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

@testset "CorrelationLengthGap: ξ = v/Δ" begin
    for (v, Δ) in ((1.0, 0.5), (2.3, 0.1), (0.7, 1.4))
        @test check(CorrelationLengthGap(); ξ=v / Δ, v=v, Δ=Δ, atol=1e-13)
        # affine in ξ and v — those solve; Δ (1/Δ) does not
        @test solve(CorrelationLengthGap(), Val(:ξ); v=v, Δ=Δ) ≈ v / Δ
        @test solve(CorrelationLengthGap(), Val(:v); ξ=v / Δ, Δ=Δ) ≈ v
    end
    @test !check(CorrelationLengthGap(); ξ=1.0, v=1.0, Δ=1.0 + 1e-3)
end

@testset "generic solve refuses a non-affine variable that blows up at a probe" begin
    # ξ = v/Δ is 1/Δ: solving for Δ would probe Δ = 0 (⇒ Inf); the guard
    # must REFUSE, not silently return NaN.
    err = try
        solve(CorrelationLengthGap(), Val(:Δ); ξ=2.0, v=1.0)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("not affine", err.msg)
end

@testset "DynamicalScaling: d(lnΔ)/d(lnξ) = −z reads off the dynamical exponent" begin
    for z in (1.0, 2.0, 0.5)
        A = 3.0
        Δ(ξ) = A * ξ^(-z)                      # synthetic gap–ξ power law
        ξ0, h = 10.0, 1e-5
        slope = (log(Δ(ξ0 * exp(h))) - log(Δ(ξ0 * exp(-h)))) / (2h)   # d lnΔ/d lnξ
        @test check(DynamicalScaling(); dlogΔ_dlogξ=slope, z=z, atol=1e-6)
        @test solve(DynamicalScaling(), Val(:z); dlogΔ_dlogξ=slope) ≈ z atol = 1e-6
    end
    # z = 1 is the Lorentz-invariant (relativistic) value
    @test check(DynamicalScaling(); dlogΔ_dlogξ=-1.0, z=1.0)
end
