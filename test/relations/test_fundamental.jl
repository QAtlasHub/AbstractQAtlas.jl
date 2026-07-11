# Fundamental thermodynamic equations vs exact closed forms.
#
# Two-level system (E ∈ {0, ε}): every potential has a closed form, so
# the algebraic relations are checked exactly and the derivative-form
# relations are checked against INDEPENDENT finite-difference routes
# (derivative vs algebra — two genuinely different computations).

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

# closed forms (total quantities, N = 1)
f2_Z(β, ε) = 1 + exp(-β * ε)
f2_F(β, ε) = -log(f2_Z(β, ε)) / β
f2_U(β, ε) = ε * exp(-β * ε) / f2_Z(β, ε)
f2_S(β, ε) = β * (f2_U(β, ε) - f2_F(β, ε))    # entropy via Legendre (exact)

@testset "FreeEnergyFromZ: statistical definition" begin
    ε = 1.3
    for T in (0.4, 1.0, 3.0)
        β = 1 / T
        F = f2_F(β, ε)
        @test check(FreeEnergyFromZ(); f=F, Z=f2_Z(β, ε), β=β, atol=1e-14)
        @test solve(FreeEnergyFromZ(), Val(:f); Z=f2_Z(β, ε), β=β) ≈ F atol = 1e-14
        # Z back-solve round-trips through the log/exp pair
        @test solve(FreeEnergyFromZ(), Val(:Z); f=F, β=β) ≈ f2_Z(β, ε) rtol = 1e-12
        # T-form kwarg equals β-form
        @test solve(FreeEnergyFromZ(), Val(:f); Z=f2_Z(β, ε), T=T) ≈ F atol = 1e-14
    end
end

@testset "FreeEnergyFromZ: per-site granularity (N)" begin
    # two independent, identical two-level systems: Z_total = Z², and the
    # per-site free energy must equal the single-system total.
    ε = 0.7
    β = 1.25
    Z1 = f2_Z(β, ε)
    f_per_site = solve(FreeEnergyFromZ(), Val(:f); Z=Z1^2, β=β, N=2)
    @test f_per_site ≈ f2_F(β, ε) atol = 1e-14
end

@testset "FreeEnergyLegendre: F = U − TS (exact algebra)" begin
    ε = 1.3
    for T in (0.4, 1.0, 3.0)
        β = 1 / T
        F, U, S = f2_F(β, ε), f2_U(β, ε), f2_S(β, ε)
        @test check(FreeEnergyLegendre(); F=F, U=U, S=S, β=β, atol=1e-14)
        @test solve(FreeEnergyLegendre(), Val(:F); U=U, S=S, β=β) ≈ F atol = 1e-14
        @test solve(FreeEnergyLegendre(), Val(:U); F=F, S=S, β=β) ≈ U atol = 1e-14
        @test solve(FreeEnergyLegendre(), Val(:S); F=F, U=U, T=T) ≈ S atol = 1e-14
        # exact-arithmetic contract: rational inputs stay rational
        r = residual(FreeEnergyLegendre(); F=1//2, U=3//2, S=1//1, β=1//1)
        @test r isa Rational
        @test r == 0//1
    end
end

@testset "EntropyResponse: S == −dF/dT (independent derivative route)" begin
    ε = 1.3
    for T in (0.5, 1.0, 2.0)
        β = 1 / T
        h = 1e-5 * T
        dF_dT = (f2_F(1 / (T + h), ε) - f2_F(1 / (T - h), ε)) / (2h)
        # derivative route vs the Legendre-algebra route for S
        @test check(EntropyResponse(); S=f2_S(β, ε), dF_dT=dF_dT, atol=1e-7)
        @test solve(EntropyResponse(), Val(:S); dF_dT=dF_dT) ≈ f2_S(β, ε) atol = 1e-7
    end
end

@testset "GibbsHelmholtz: U == ∂(βF)/∂β (independent derivative route)" begin
    ε = 1.3
    for β in (0.3, 1.0, 2.5)
        h = 1e-6 * β
        dβF = (((β + h) * f2_F(β + h, ε)) - ((β - h) * f2_F(β - h, ε))) / (2h)
        @test check(GibbsHelmholtz(); U=f2_U(β, ε), dβF_dβ=dβF, atol=1e-7)
        @test solve(GibbsHelmholtz(), Val(:U); dβF_dβ=dβF) ≈ f2_U(β, ε) atol = 1e-7
    end
end

@testset "the web closes: Z → F → (S, U) → F reconciles three routes" begin
    # Z-definition, Legendre algebra, and both derivative routes must
    # reconcile on the same state point — the QAtlas identities-plane
    # reconciliation expressed through this package's relations.
    ε = 0.9
    T = 1.7
    β = 1 / T
    F = solve(FreeEnergyFromZ(), Val(:f); Z=f2_Z(β, ε), β=β)
    h = 1e-5
    S = solve(
        EntropyResponse(),
        Val(:S);
        dF_dT=(f2_F(1 / (T + h), ε) - f2_F(1 / (T - h), ε)) / (2h),
    )
    U = solve(
        GibbsHelmholtz(),
        Val(:U);
        dβF_dβ=(((β + h) * f2_F(β + h, ε)) - ((β - h) * f2_F(β - h, ε))) / (2h),
    )
    @test check(FreeEnergyLegendre(); F=F, U=U, S=S, β=β, atol=1e-6)
end
