# Dynamical / spectral pointwise identities vs independent constructions.
#
# Each identity is checked against a value built a DIFFERENT way than the
# identity states: Dyson against a G assembled from G₀ and Σ; the
# spectral representation against an explicit Lorentzian; detailed
# balance against an exponential; the NMR exponent against exact
# rationals.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

@testset "Dyson equation (complex, pointwise)" begin
    # build the FULL G from bare G₀ and Σ, then the identity must close
    for (G0, Σ) in [(1 / (0.5 + 0.1im), 0.2 - 0.3im), (1 / (-1.2 + 0.4im), -0.5 + 0.2im)]
        G = 1 / (1 / G0 - Σ)
        @test abs(residual(Dyson(); G=G, G0=G0, Σ=Σ)) < 1e-14
        # solve for the self-energy from G and G₀ (generic affine solve)
        @test solve(Dyson(), Val(:Σ); G=G, G0=G0) ≈ Σ
    end
end

@testset "spectral representation A = −Im G^R/π (Lorentzian)" begin
    # a single pole at ε with width γ: G^R(ω) = 1/(ω − ε + iγ),
    # A(ω) = (1/π) γ/((ω−ε)² + γ²) — an independent construction.
    ε, γ = 0.3, 0.05
    for ω in (0.0, 0.3, 0.7)
        GR = 1 / (ω - ε + im * γ)
        A = (1 / π) * γ / ((ω - ε)^2 + γ^2)
        @test check(SpectralFromGreens(); A=A, G=GR, atol=1e-13)
    end
    # and the sum rule: ∫ A dω = 1 for the Lorentzian (analytic)
    @test check(SpectralSumRule(); spectral_integral=1.0, atol=0)
end

@testset "detailed balance S(q,−ω) = e^{−βω} S(q,ω)" begin
    for T in (0.5, 2.0), ω in (0.4, 1.3)
        β = 1 / T
        Splus = 3.7                      # arbitrary S(q, ω)
        Sminus = exp(-β * ω) * Splus     # the balanced partner
        @test check(DetailedBalance(); S_plus=Splus, S_minus=Sminus, ω=ω, β=β, atol=1e-13)
        # T-form normalizes identically
        @test check(DetailedBalance(); S_plus=Splus, S_minus=Sminus, ω=ω, T=T, atol=1e-13)
        # solve for the anti-Stokes side
        @test solve(DetailedBalance(), Val(:S_minus); S_plus=Splus, ω=ω, β=β) ≈ Sminus
    end
end

@testset "NMR exponent θ_NMR = 2Δ_op − 1 (exact rationals)" begin
    # 1D transverse-field Ising QCP: Δ_σ = 1/8 ⟹ θ_NMR = −3/4, exactly
    @test residual(NMRExponent(); θ_NMR=-3 // 4, Δ_op=1 // 8) == 0 // 1
    @test solve(NMRExponent(), Val(:θ_NMR); Δ_op=1 // 8) == -3 // 4
    @test solve(NMRExponent(), Val(:Δ_op); θ_NMR=-3 // 4) == 1 // 8
    @test solve(NMRExponent(), Val(:θ_NMR); Δ_op=1 // 8) isa Rational   # exactness
    @test !check(NMRExponent(); θ_NMR=0 // 1, Δ_op=1 // 8)
end

@testset "structure-factor sum rule: S(q) = ∫ S(q,ω) dω/(2π)" begin
    @test residual(StaticFromDynamicalStructureFactor(); Sq=1.5, sqw_integral=1.5) == 0
    @test check(StaticFromDynamicalStructureFactor(); Sq=2.0, sqw_integral=2.0)
    @test !check(StaticFromDynamicalStructureFactor(); Sq=2.0, sqw_integral=1.0)
    @test solve(StaticFromDynamicalStructureFactor(), Val(:Sq); sqw_integral=0.7) ≈ 0.7
end

@testset "dynamical FDT reproduces (convention-free) detailed balance" begin
    # build S(±ω) from χ''(±ω) via the FDT, then the INDEPENDENT detailed-
    # balance relation must hold — a cross-check between two relations.
    for T in (0.5, 2.0), ω in (0.4, 1.3), χpp in (0.6, 2.5)
        β = 1 / T
        Sp = solve(DynamicalFDT(), Val(:S); χpp=χpp, ω=ω, β=β)
        Sm = solve(DynamicalFDT(), Val(:S); χpp=(-χpp), ω=(-ω), β=β)   # χ'' is odd
        @test isapprox(Sm / Sp, exp(-β * ω); rtol=1e-12)
        @test check(DetailedBalance(); S_plus=Sp, S_minus=Sm, ω=ω, β=β, atol=1e-12)
        # T-form kwarg agrees with β-form
        @test solve(DynamicalFDT(), Val(:S); χpp=χpp, ω=ω, T=T) ≈ Sp
    end
end

@testset "Kramers–Kronig: exact arithmetic + convention vs a Hilbert transform" begin
    # exactness of the relation itself (the 1/π factor and the signs)
    for (Reχ, Imχ) in ((0.7, -0.4), (-1.3, 2.1))
        @test check(KramersKronigReal(); Reχ=Reχ, pv_imag=π * Reχ, atol=1e-12)
        @test check(KramersKronigImag(); Imχ=Imχ, pv_real=(-π * Imχ), atol=1e-12)
        @test solve(KramersKronigReal(), Val(:Reχ); pv_imag=π * Reχ) ≈ Reχ
        @test solve(KramersKronigImag(), Val(:Imχ); pv_real=(-π * Imχ)) ≈ Imχ
    end

    # INDEPENDENT convention check: the retarded Green's function of a single
    # level (ε = 0, width Γ) is an exact Kramers–Kronig pair; feed the
    # numerically Hilbert-transformed part and confirm the OTHER part.
    Γ = 1.0
    ReG(ω) = ω / (ω^2 + Γ^2)
    ImG(ω) = -Γ / (ω^2 + Γ^2)
    # P∫ g(ω')/(ω'−ω) dω' by pole subtraction:
    # ∫[g(ω')−g(ω)]/(ω'−ω) dω' + g(ω)·ln|(L−ω)/(L+ω)| (bracket smooth at ω'=ω)
    function hilbert_pv(g, ω; L=5000.0, n=2_000_000)
        h = 2L / n
        s = 0.0
        gω = g(ω)
        for k in 0:n
            ωp = -L + k * h
            s +=
                ((k == 0 || k == n) ? 0.5 : 1.0) * (ωp == ω ? 0.0 : (g(ωp) - gω) / (ωp - ω))
        end
        return s * h + gω * log((L - ω) / (L + ω))
    end
    for ω in (0.3, 1.2, -0.7)
        # χ'(ω) = (1/π) P∫ χ''/(ω'−ω) reproduces ReG; χ''(ω) = −(1/π) P∫ χ'/(ω'−ω) reproduces ImG
        @test check(KramersKronigReal(); Reχ=ReG(ω), pv_imag=hilbert_pv(ImG, ω), atol=3e-4)
        @test check(KramersKronigImag(); Imχ=ImG(ω), pv_real=hilbert_pv(ReG, ω), atol=3e-4)
    end
end

@testset "one-call sweep picks up the spectral identities by variable name" begin
    G0 = 1 / (0.5 + 0.1im)
    Σ = 0.2 - 0.3im
    G = 1 / (1 / G0 - Σ)
    rep = relation_report((; G=G, G0=G0, Σ=Σ); atol=1e-12, domain=:spectral)
    @test length(rep) == 1
    @test rep[1].relation isa Dyson
    @test rep[1].pass
end

@testset "StaticStructureFactorFromCorrelation: S(q→0) = ∫G(r)dr" begin
    R = StaticStructureFactorFromCorrelation
    # exact residual: zero iff S(q=0) equals the supplied spatial integral
    @test residual(R(); Sq0=3 // 2, integral_G=3 // 2) == 0 // 1
    @test check(R(); Sq0=2.0, integral_G=2.0)
    @test !check(R(); Sq0=2.0, integral_G=1.0)
    # affine ⇒ solve either variable
    @test solve(R(), Val(:Sq0); integral_G=5 // 1) == 5 // 1
    @test solve(R(), Val(:integral_G); Sq0=7 // 2) == 7 // 2
    # type-keyed: Sq0 is a StaticStructureFactor slot ⇒ auto-links to that quantity
    @test StaticStructureFactor in quantities(R())
    @test R() in relations_constraining(StaticStructureFactor)
    # the collision-proof bag front door
    @test residual(R(), bag(StaticStructureFactor => 2.0); integral_G=2.0) == 0.0
end

@testset "FSumRule: ∫ω S(q,ω)dω = N q²/(2m) (f-sum rule)" begin
    F = FSumRule
    # per-particle (N=1): first moment = q²/(2m); exact (Rational in ⇒ Rational out)
    @test residual(F(); first_moment=2 // 1, q=2 // 1, m=1 // 1) == 0 // 1   # 2²/(2·1) = 2
    @test check(F(); first_moment=2.0, q=2.0, m=1.0)
    @test !check(F(); first_moment=1.0, q=2.0, m=1.0)
    # N particles scale the moment
    @test residual(F(); first_moment=6 // 1, q=2 // 1, m=1 // 1, N=3) == 0 // 1   # 3·4/2 = 6
    # affine ⇒ solve the supplied moment
    @test solve(F(), Val(:first_moment); q=2 // 1, m=1 // 1) == 2 // 1
    # hand-linked to the dynamical structure factor (constrained via the supplied moment)
    @test DynamicalStructureFactor in quantities(F())
    @test F() in relations_constraining(DynamicalStructureFactor)
end
