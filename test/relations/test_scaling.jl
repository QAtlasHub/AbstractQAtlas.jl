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

# ─────────────────────────────────────────────────────────────────────────────
# The infinite-randomness exponent network, checked against the ONE exact
# exponent set that pins it: the 1D random transverse-field Ising chain,
# Iglói–Monthus Table 1, §4.1.2 (ν = 2, ψ = 1/2, ν_typ = 1, x_m = (3−√5)/4).  Each
# relation below is stated in the review for general `d`; the table is 1D.  So
# these are not restatements — each is the general law evaluated at a point it
# did not come from.
# ─────────────────────────────────────────────────────────────────────────────

const RTFIC_ν = 2//1
const RTFIC_ψ = 1//2
const RTFIC_ν_typ = 1 // 1               # Table 1 and Eq. (4.10), obtained separately
const RTFIC_x_m = (3 - sqrt(5)) / 4

@testset "the typical correlation length is the table's ν_typ, not ν" begin
    # ν(1−ψ) = 2·(1/2) = 1 exactly, and Eq. (4.10) reaches ν_typ = 1 from the
    # finite-size dependence instead — two routes, one number.
    @test residual(TypicalCorrelationLength(); ν_typ=RTFIC_ν_typ, ν=RTFIC_ν, ψ=RTFIC_ψ) ==
        0//1
    @test solve(TypicalCorrelationLength(), Val(:ν_typ); ν=RTFIC_ν, ψ=RTFIC_ψ) ==
        RTFIC_ν_typ
    @test solve(TypicalCorrelationLength(), Val(:ψ); ν_typ=RTFIC_ν_typ, ν=RTFIC_ν) ==
        RTFIC_ψ

    # ν_typ < ν whenever ψ > 0, and they coincide only at ψ = 0 — i.e. only
    # where the fixed point is not infinite-randomness.
    @test solve(TypicalCorrelationLength(), Val(:ν_typ); ν=RTFIC_ν, ψ=0//1) == RTFIC_ν
    @test solve(TypicalCorrelationLength(), Val(:ν_typ); ν=RTFIC_ν, ψ=RTFIC_ψ) < RTFIC_ν

    # A fixture that CAN disagree: ν_typ = ν is the clean-system answer and is
    # wrong here.
    @test !check(TypicalCorrelationLength(); ν_typ=RTFIC_ν, ν=RTFIC_ν, ψ=RTFIC_ψ)
end

@testset "the Griffiths exponent's divergence rate reproduces the exact 1D law" begin
    # 1D has the closed-form Griffiths result 1/z = 2|δ| (stated with Eq. (4.51),
    # §4.4.2), so z = 1/(2|δ|) and d(ln z)/d(ln|δ|) = −1 identically.
    # The general law says the rate is −νψ; with the table's ν, ψ that is −1.
    # Two independent routes to the same number.
    @test solve(GriffithsExponentDivergence(), Val(:dlogz_dlogδ); ν=RTFIC_ν, ψ=RTFIC_ψ) ==
        -1//1

    # ...and measure that slope numerically off z = 1/(2|δ|) rather than
    # asserting the algebra twice.
    δs = 10.0 .^ range(-4, -2; length=25)
    zs = 1 ./ (2 .* δs)
    lx, ly = log.(δs), log.(zs)
    slope =
        ((lx .- sum(lx) / length(lx))' * (ly .- sum(ly) / length(ly))) /
        sum(abs2, lx .- sum(lx) / length(lx))
    @test check(GriffithsExponentDivergence(); dlogz_dlogδ=slope, ν=2.0, ψ=0.5, atol=1e-10)

    # The clean case is a different law, not this one with ψ = 0: at ψ = 0 the
    # rate is zero, i.e. z does not drift — which is exactly what a conventional
    # critical point does and what this relation must NOT claim for the random one.
    @test !check(GriffithsExponentDivergence(); dlogz_dlogδ=slope, ν=2.0, ψ=0.0, atol=1e-6)
end

@testset "the moment-growth exponent returns the golden mean" begin
    # Eq. (A.21) at d = 1, x_m = (3−√5)/4, ψ = 1/2 must give (1+√5)/2 — the value
    # Eq. (3.18) reaches by a different route, here from two other Table 1 entries.
    φ = solve(ActivatedMomentGrowth(), Val(:φ); d=1.0, x_m=RTFIC_x_m, ψ=0.5)
    @test φ ≈ (1 + sqrt(5)) / 2 rtol = 1e-14
    # It is not 2 and not 1.5 — pin that the check discriminates.
    @test !check(ActivatedMomentGrowth(); φ=2.0, d=1.0, x_m=RTFIC_x_m, ψ=0.5, atol=1e-6)
end

@testset "the two Griffiths observables read the same z" begin
    # χ ∼ T^{−1+d/z} and c_V ∼ T^{d/z} are different measurements; a z fitted
    # from one alone is a fit, the pair is a check.  d = 1, z = 4.
    d, z = 1//1, 4//1
    dlogχ, dlogc = -1 + d // z, d // z
    @test residual(GriffithsSusceptibility(); dlogχ_dlogT=dlogχ, d=d, z=z) == 0//1
    @test residual(GriffithsSpecificHeat(); dlogc_dlogT=dlogc, d=d, z=z) == 0//1
    @test solve(GriffithsSusceptibility(), Val(:z); dlogχ_dlogT=dlogχ, d=d) == z
    @test solve(GriffithsSpecificHeat(), Val(:z); dlogc_dlogT=dlogc, d=d) == z

    # z = d is the line inside the Griffiths phase where χ stops diverging.
    @test solve(GriffithsSusceptibility(), Val(:dlogχ_dlogT); d=1//1, z=1//1) == 0//1
    @test solve(GriffithsSusceptibility(), Val(:dlogχ_dlogT); d=1//1, z=4//1) < 0

    # The two observables disagree if fed the same slope — they must not be
    # the same relation wearing two names.
    @test !check(GriffithsSpecificHeat(); dlogc_dlogT=dlogχ, d=d, z=z)
end

@testset "the infinite-randomness relations key on the reductions" begin
    # `Typical` and `DisorderAveraged` are separate quantities so a bag cannot mix
    # them, and each relation here is about one of them. `also_constrains` is
    # family-erased by design, so the graph carries "a disorder reduction" and not
    # which observable was reduced — coarser than the key, and the same coarseness
    # `TypicalBelowAverage` already has. Lookup is `Q <: T`, so the wrapped type
    # still finds it.
    for (rel, wrapped) in (
        (GriffithsSusceptibility(), DisorderAveraged{Susceptibility}),
        (GriffithsSpecificHeat(), DisorderAveraged{SpecificHeat}),
        (TypicalCorrelationLength(), Typical{CorrelationLength}),
        (TypicalCorrelationLength(), DisorderAveraged{CorrelationLength}),
        (ActivatedDynamicalScaling(), Typical{MassGap}),
    )
        @test rel in relations_constraining(wrapped)
    end
    @test GriffithsExponentDivergence() in relations_constraining(DynamicalExponent)

    # The clean quantity must NOT reach them: a clean susceptibility has no
    # Griffiths singularity, and there is only one correlation length in a clean
    # system — which is exactly what `TypicalCorrelationLength` exists to deny.
    for (rel, clean) in (
        (GriffithsSusceptibility(), Susceptibility{(:z, :z)}),
        (GriffithsSpecificHeat(), SpecificHeat),
        (TypicalCorrelationLength(), CorrelationLength),
    )
        @test !(rel in relations_constraining(clean))
    end

    @test ActivatedExponent in quantities(TypicalCorrelationLength())
    @test ActivatedExponent in quantities(GriffithsExponentDivergence())
    @test ActivatedExponent in quantities(ActivatedMomentGrowth())
end

@testset "a plain exponent set does not silently pick up the disorder laws" begin
    # `exponents_consistent` sweeps every :scaling relation.  The five added
    # here need ψ / z / ν_typ / φ / x_m, none of which a clean (α,β,γ,δ,ν,η) set
    # carries — so the gate must be unchanged for the sets it already served.
    @test exponents_consistent(ISING2D; d=2)
    names = keys(exponent_residuals(ISING2D; d=2))
    for n in (:typicalcorrelationlength, :griffithssusceptibility, :activatedmomentgrowth)
        @test !(n in names)
    end
end

@testset "the finite-size gap tells a CFT point from an infinite-randomness one" begin
    localslope(y, x, k) = (y[k + 1] - y[k]) / (x[k + 1] - x[k])
    Ls = (16.0, 32.0, 64.0, 128.0, 256.0)
    lnL = log.(collect(Ls))

    # Two gaps over the SAME sizes: Cardy's 2πvx/L on a periodic chain, and the
    # activated exp(-c L^ψ) of an open random one.
    v, x, c, ψtrue = 1.0, 1 // 8, 0.7, 0.5
    Δ_cft = [2π * v * x / L for L in Ls]
    Δ_irfp = [exp(-c * L^ψtrue) for L in Ls]

    # Each law is exact on its own data.
    @test all(L -> check(FiniteSizeGap(); gap=2π * v * x / L, x=x, v=v, L=L), Ls)
    ψ_irfp = [localslope(log.(-log.(Δ_irfp)), lnL, k) for k in 1:4]
    @test all(s -> isapprox(s, ψtrue; atol=1e-12), ψ_irfp)
    @test check(
        ActivatedFiniteSizeScaling(); dloglogO_dlogL=ψ_irfp[end], ψ=ψtrue, atol=1e-12
    )
    @test solve(ActivatedFiniteSizeScaling(), Val(:ψ); dloglogO_dlogL=ψ_irfp[end]) ≈ ψtrue

    # And each rejects the other's, which is what makes one sweep in L a test.
    # Reading a scaling dimension off the activated gap gives a different answer
    # at every size, and the drift is derived rather than eyeballed:
    # x(L) = exp(-c L^ψ)·L/(2πv), so the ends of the sweep differ by exactly
    # exp(c(L_max^ψ - L_min^ψ))·L_min/L_max, which diverges as the sweep grows.
    x_from_irfp = [Δ * L / (2π * v) for (Δ, L) in zip(Δ_irfp, Ls)]
    @test !check(FiniteSizeGap(); gap=Δ_irfp[end], x=x_from_irfp[1], v=v, L=Ls[end])
    @test x_from_irfp[1] / x_from_irfp[end] ≈
        exp(c * (Ls[end]^ψtrue - Ls[1]^ψtrue)) * Ls[1] / Ls[end]
    @test issorted(x_from_irfp; rev=true)
    # Conversely the CFT gap has no constant ψ: its apparent one decays as
    # 1/ln L rather than sitting still.
    ψ_cft = [localslope(log.(-log.(Δ_cft)), lnL, k) for k in 1:4]
    @test !check(ActivatedFiniteSizeScaling(); dloglogO_dlogL=ψ_cft[end], ψ=ψtrue, atol=0.1)
    @test issorted(ψ_cft; rev=true)

    @test variable_types(ActivatedFiniteSizeScaling()) == (ActivatedExponent,)
    @test ActivatedFiniteSizeScaling() in relations_constraining(ActivatedExponent)
end
