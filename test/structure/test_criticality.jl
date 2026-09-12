# The quantity ⇄ exponent correspondence, and the forms DERIVED from it.
#
# Independent expectations: the correspondence must reproduce the exact
# textbook power laws (M∼|t|^β, χ∼|t|^{-γ}, ξ∼|t|^{-ν}, C∼|t|^{-α}) and
# the FSS combinations (χ_max∼L^{γ/ν}, M∼L^{-β/ν}, ξ∼L) — and it must do
# so by LOOKING UP the exponent, not by the caller passing it.

using AbstractQAtlas
using AbstractQAtlas: critical_scaling, singular_form, fss_size_exponent, fss_peak

const EXPS = (α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4)   # 2D Ising, exact

@testset "correspondence: which exponent governs which quantity" begin
    @test critical_scaling(SpontaneousMagnetization()) == CriticalScaling(:β, +1)
    @test critical_scaling(Susceptibility(:z, :z)) == CriticalScaling(:γ, -1)
    @test critical_scaling(Susceptibility(:x, :x)) == CriticalScaling(:γ, -1)  # any axis
    @test critical_scaling(SpecificHeat()) == CriticalScaling(:α, -1)
    @test critical_scaling(CorrelationLength()) == CriticalScaling(:ν, -1)
    # quantities with no reduced-temperature critical law
    @test critical_scaling(PartitionFunction()) === nothing
    @test critical_scaling(FreeEnergy()) === nothing
    # the field-driven / distance-driven laws have their own accessors
    @test critical_isotherm(SpontaneousMagnetization()) == :δ
    @test correlation_decay(SpinCorrelation(:z, :z)) == :η
end

@testset "singular_form derived from the correspondence (exact rationals)" begin
    # exponent value + sign both come from the correspondence, not the call
    @test singular_form(SpontaneousMagnetization(), -1 // 1; exponents=EXPS) == 1  # |−1|^{1/8}
    # M ∼ |t|^{+β}: order parameter grows as |t| grows on the ordered side
    m1 = singular_form(SpontaneousMagnetization(), -0.01; exponents=EXPS)
    m2 = singular_form(SpontaneousMagnetization(), -0.04; exponents=EXPS)
    @test m2 > m1
    @test m1 ≈ 0.01^(1 / 8)
    # χ ∼ |t|^{-γ}: diverges as t → 0
    x1 = singular_form(Susceptibility(:z, :z), 0.01; exponents=EXPS)
    x2 = singular_form(Susceptibility(:z, :z), 0.005; exponents=EXPS)
    @test x2 > x1
    @test x1 ≈ 0.01^(-7 / 4)
    # ξ ∼ |t|^{-ν}, C ∼ |t|^{-α}
    @test singular_form(CorrelationLength(), 0.02; exponents=EXPS) ≈ 0.02^(-1.0)
    @test singular_form(SpecificHeat(), 0.02; exponents=EXPS) == 1  # α = 0 ⇒ |t|^0
    # refused for a non-critical quantity
    @test_throws ErrorException singular_form(PartitionFunction(), 0.1; exponents=EXPS)
end

@testset "FSS size exponent derived (the γ/ν that used to live in a comment)" begin
    # χ_max ∼ L^{+γ/ν}
    @test fss_size_exponent(Susceptibility(:z, :z); exponents=EXPS) == 7 // 4
    # M(Tc) ∼ L^{-β/ν}
    @test fss_size_exponent(SpontaneousMagnetization(); exponents=EXPS) == -1 // 8
    # ξ ∼ L^{1} — correlation length saturates at the system size
    @test fss_size_exponent(CorrelationLength(); exponents=EXPS) == 1 // 1
    # C_max ∼ L^{α/ν}
    @test fss_size_exponent(SpecificHeat(); exponents=EXPS) == 0 // 1
    # exactness preserved through the derivation (Rational in ⇒ Rational out)
    @test fss_size_exponent(Susceptibility(:z, :z); exponents=EXPS) isa Rational
end

@testset "fss_peak round trip: synthetic data recovers the derived exponent" begin
    ratio = fss_size_exponent(Susceptibility(:z, :z); exponents=EXPS)   # 7/4, derived
    Ls = [8, 16, 32, 64]
    χmax = [2.31 * fss_peak(Susceptibility(:z, :z), L; exponents=EXPS) for L in Ls]
    # least-squares log-log slope must recover the *derived* ratio
    mx, my = sum(log.(Ls)) / 4, sum(log.(χmax)) / 4
    slope = sum((log.(Ls) .- mx) .* (log.(χmax) .- my)) / sum((log.(Ls) .- mx) .^ 2)
    @test isapprox(slope, float(ratio); atol=1e-12)
end

@testset "collapse_coordinates: quantity-driven, exponents from the atlas" begin
    Tc = 2.269185314213022
    # pivot: x = 0 at T = Tc for every L; scale = L^{ρ} with ρ = −fss_size_exponent
    for L in (8, 16, 32)
        c = collapse_coordinates(SpontaneousMagnetization(), Tc, L, Tc; exponents=EXPS)
        @test c.x == 0.0
        @test c.scale ≈ float(L)^(1 / 8)          # ρ = +β/ν = 1/8
    end
    cχ = collapse_coordinates(Susceptibility(:z, :z), Tc + 0.01, 16, Tc; exponents=EXPS)
    @test cχ.x ≈ 0.01 * 16.0                        # (T−Tc)·L^{1/ν}
    @test cχ.scale ≈ 16.0^(-7 / 4)                  # ρ = −γ/ν
    # perfect two-size collapse of a synthetic observable
    f(x) = exp(-x^2)
    ν = 1.0
    ρ = -float(fss_size_exponent(SpontaneousMagnetization(); exponents=EXPS))  # 1/8
    for L1 in (8, 64), L2 in (16, 32)
        xs = 0.37
        T1 = Tc + xs * L1^(-1 / ν)
        T2 = Tc + xs * L2^(-1 / ν)
        O1 =
            L1^(-ρ) * f(
                collapse_coordinates(
                    SpontaneousMagnetization(), T1, L1, Tc; exponents=EXPS
                ).x,
            )
        O2 =
            L2^(-ρ) * f(
                collapse_coordinates(
                    SpontaneousMagnetization(), T2, L2, Tc; exponents=EXPS
                ).x,
            )
        s1 =
            collapse_coordinates(SpontaneousMagnetization(), T1, L1, Tc; exponents=EXPS).scale
        s2 =
            collapse_coordinates(SpontaneousMagnetization(), T2, L2, Tc; exponents=EXPS).scale
        @test isapprox(O1 * s1, O2 * s2; rtol=1e-12)
    end
end

# Igloi-Monthus Table 1 (§4.1.2): the random transverse-field Ising chain.
const SDRG = (β=(3 - sqrt(5)) / 2, β_s=1 // 1, ν=2 // 1, ψ=1 // 2)

@testset "at an infinite-randomness point a bare quantity names no size law" begin
    m = SurfaceMagnetization()
    @test critical_scaling(m) == CriticalScaling(:β_s, +1)

    # The disorder average is an ordinary power of L, and the power it returns
    # is the review's own surface dimension: -x_m_s = -1/2, Eq. (4.7).
    @test fss_size_exponent(DisorderAveraged(m); exponents=SDRG) == -1 // 2
    # Same for the bulk order parameter: -x_m = -(3-√5)/4, Table 1.
    @test fss_size_exponent(DisorderAveraged(SpontaneousMagnetization()); exponents=SDRG) ≈
        -(3 - sqrt(5)) / 4

    # The typical value has no power of L at all, and the refusal says that
    # rather than "this quantity has no entry".
    @test_throws "no power of L" fss_size_exponent(Typical(m); exponents=SDRG)
    @test_throws "ActivatedFiniteSizeScaling" fss_size_exponent(Typical(m); exponents=SDRG)
    @test_throws "ActivatedFiniteSizeScaling" fss_peak(Typical(m), 64; exponents=SDRG)
    @test_throws "ActivatedFiniteSizeScaling" collapse_coordinates(
        Typical(m), 0.1, 64, 0.0; exponents=SDRG
    )

    # A bare quantity is the trap: the answer -1/2 is right for the average and
    # wrong for the typical, so once ψ is present the caller has to say which.
    @test_throws "DisorderAveraged" fss_size_exponent(m; exponents=SDRG)
    # Nothing is refused where there is nothing to disambiguate: drop ψ, or set
    # it to zero (not an infinite-randomness point), and the same call answers.
    @test fss_size_exponent(m; exponents=(β_s=1 // 1, ν=2 // 1)) == -1 // 2
    @test fss_size_exponent(m; exponents=(β_s=1 // 1, ν=2 // 1, ψ=0 // 1)) == -1 // 2
end
