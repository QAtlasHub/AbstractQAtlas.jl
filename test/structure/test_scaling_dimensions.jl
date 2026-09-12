# The RG-eigenvalue origin of the critical exponents: `critical_exponents`
# derives (α,β,γ,δ,ν,η) from (y_t, y_h, d), and the four scaling laws in
# relations/scaling.jl become CONSEQUENCES — their residual is identically 0
# for the derived set, at any eigenvalues.  The test checks the derivation
# against independent expectations (the scaling @relations themselves, known
# exponent tables, the inverse map), NOT against its own output.

using AbstractQAtlas
using AbstractQAtlas:
    ScalingDimensions,
    InfiniteRandomness,
    critical_exponents,
    critical_exponent,
    scaling_dimensions,
    infinite_randomness,
    residual,
    Rushbrooke,
    Widom,
    Fisher,
    Josephson,
    exponents_consistent,
    exponent_residuals

# a spread of EXACT rational RG data (various y_t, y_h, d), including the two
# physical fixed points 2D-Ising (1, 15//8, 2) and 3D-percolation-like sets
const _RG_SETS = (
    (1 // 1, 15 // 8, 2),      # 2D Ising (exact)
    (3 // 2, 7 // 4, 3),
    (2 // 1, 9 // 5, 3),
    (5 // 4, 11 // 6, 2),
    (4 // 3, 21 // 11, 3),
)

@testset "the four scaling laws are IDENTITIES in (y_t, y_h, d)" begin
    # feed the DERIVED exponents to the independent scaling @relations — every
    # residual must be exactly 0//1, for every eigenvalue set.  This is the
    # structural claim: Rushbrooke/Widom/Fisher/Josephson are not axioms but
    # consequences of two-eigenvalue homogeneity + hyperscaling.
    for (yt, yh, dd) in _RG_SETS
        e = critical_exponents(ScalingDimensions(yt, yh, dd))
        @test residual(Rushbrooke(); α=e.α, β=e.β, γ=e.γ) == 0 // 1
        @test residual(Widom(); β=e.β, γ=e.γ, δ=e.δ) == 0 // 1
        @test residual(Fisher(); γ=e.γ, ν=e.ν, η=e.η) == 0 // 1
        @test residual(Josephson(); α=e.α, ν=e.ν, d=dd) == 0 // 1
        # and the domain-wide gate agrees
        @test exponents_consistent(e; d=dd)
        # exactness: Rational eigenvalues ⇒ Rational exponents
        @test all(v -> v isa Rational, values(e))
    end
end

@testset "2D Ising exponents fall out of (1, 15//8, 2)" begin
    e = critical_exponents(ScalingDimensions(1 // 1, 15 // 8, 2))
    @test e == (α=0 // 1, β=1 // 8, γ=7 // 4, δ=15 // 1, ν=1 // 1, η=1 // 4)
    # single-exponent accessor agrees with the full set
    for name in (:α, :β, :γ, :δ, :ν, :η)
        @test critical_exponent(name, ScalingDimensions(1 // 1, 15 // 8, 2)) == e[name]
    end
end

@testset "mean-field / Gaussian eigenvalues at the upper critical dimension" begin
    # Gaussian fixed point: y_t = 2, y_h = 1 + d/2.  At the Ising upper critical
    # dimension d = 4 this gives the classical (Landau) exponents AND — since
    # hyperscaling holds exactly at d_upper — still satisfies Josephson.
    s = ScalingDimensions(2 // 1, 1 + 4 // 2, 4)   # y_h = 3, d = 4
    e = critical_exponents(s)
    @test e == (α=0 // 1, β=1 // 2, γ=1 // 1, δ=3 // 1, ν=1 // 2, η=0 // 1)
    @test exponents_consistent(e; d=4)
end

@testset "inverse map: (ν, η, d) reconstructs the eigenvalues and the rest" begin
    for (yt, yh, dd) in _RG_SETS
        e = critical_exponents(ScalingDimensions(yt, yh, dd))
        # recover the eigenvalues from just ν and η at dimension d …
        s = scaling_dimensions(; ν=e.ν, η=e.η, d=dd)
        @test s.y_t == yt
        @test s.y_h == yh
        # … and the FULL exponent set round-trips (δ, β, γ, α all reconstructed
        # from ν, η, d alone — the two-eigenvalue structure)
        @test critical_exponents(s) == e
    end
end

@testset "float eigenvalues propagate (numerical fixed point)" begin
    s = ScalingDimensions(1.0, 1.875, 2.0)
    e = critical_exponents(s)
    @test e.β ≈ 0.125
    @test e.γ ≈ 1.75
    @test e.δ ≈ 15.0
    @test e.η ≈ 0.25
    # laws still hold to floating tolerance
    @test residual(Rushbrooke(); α=e.α, β=e.β, γ=e.γ) ≈ 0 atol = 1e-12
    @test residual(Fisher(); γ=e.γ, ν=e.ν, η=e.η) ≈ 0 atol = 1e-12
end

# An infinite-randomness fixed point has the same DECLARE-ONCE structure with a
# different independent set: (ψ, ν, x_m) plus d. The three relations it feeds
# are stated separately in relations/scaling.jl, so they are the independent
# expectation here, exactly as the four clean laws are above.
const _IRFP_SETS = (
    (1 // 2, 2 // 1, 1 // 4, 1),      # RTFIC ψ and ν, rational x_m stand-in
    (1 // 2, 2 // 1, 1 // 3, 2),
    (1 // 3, 3 // 2, 1 // 5, 2),
    (2 // 3, 5 // 2, 2 // 7, 3),
    (1 // 4, 4 // 1, 1 // 8, 1),
)

@testset "the infinite-randomness relations are IDENTITIES in (ψ, ν, x_m, d)" begin
    for (ψ, ν, x_m, dd) in _IRFP_SETS
        e = critical_exponents(InfiniteRandomness(ψ, ν, x_m, dd))
        @test residual(OrderParameterDimension(); β=e.β, ν=e.ν, x_m=e.x_m) == 0 // 1
        @test residual(TypicalCorrelationLength(); ν_typ=e.ν_typ, ν=e.ν, ψ=e.ψ) == 0 // 1
        @test residual(ActivatedMomentGrowth(); φ=e.φ, d=dd, x_m=e.x_m, ψ=e.ψ) == 0 // 1
        # Rational in ⇒ Rational out, as for the clean set
        @test all(v -> v isa Rational, values(e))
        # and the gate agrees without being told d a second time
        @test exponents_consistent(InfiniteRandomness(ψ, ν, x_m, dd))
    end
end

@testset "RTFIC Table 1 falls out of (1//2, 2, x_m, d = 1)" begin
    # Igloi-Monthus Table 1 (§4.1.2). x_m is irrational, so ν_typ and the
    # rational-valued entries are exact and β, φ are compared at machine
    # precision against the table's own closed forms.
    x_m = (3 - sqrt(5)) / 4
    e = critical_exponents(InfiniteRandomness(1 // 2, 2 // 1, x_m, 1))
    @test e.ν_typ == 1.0                        # Eq. (4.10)
    @test e.β ≈ (3 - sqrt(5)) / 2               # Table 1, β = ν x_m
    @test e.φ ≈ (1 + sqrt(5)) / 2               # the golden mean, Eq. (A.21)
    for name in (:β, :ν, :ν_typ, :ψ, :x_m, :φ)
        @test critical_exponent(name, InfiniteRandomness(1 // 2, 2 // 1, x_m, 1)) == e[name]
    end
end

@testset "`d` in the struct is what makes the Euclidean mix-up one decision" begin
    # The trap: an atlas hands out a `d` kwarg meaning the EUCLIDEAN dimension
    # (2 for a 1D quantum chain), while every infinite-randomness relation
    # reads the SPATIAL one. The wrong answer does not look wrong.
    x_m = (3 - sqrt(5)) / 4
    spatial = critical_exponents(InfiniteRandomness(1 // 2, 2 // 1, x_m, 1))
    euclidean = critical_exponents(InfiniteRandomness(1 // 2, 2 // 1, x_m, 2))
    @test spatial.φ ≈ (1 + sqrt(5)) / 2
    @test euclidean.φ ≈ (1 + sqrt(5)) / 2 + 2   # 3.618, and nothing flags it
    # Only φ moves: d enters no other derived exponent, which is exactly why a
    # wrong d survives a partial check.
    for name in (:β, :ν, :ν_typ, :ψ, :x_m)
        @test spatial[name] == euclidean[name]
    end
    # The struct-aware gate has no second place to state d, so both are internally
    # consistent; what the struct buys is that the choice is made ONCE.
    @test exponents_consistent(InfiniteRandomness(1 // 2, 2 // 1, x_m, 1); atol=1e-15)
    @test exponents_consistent(InfiniteRandomness(1 // 2, 2 // 1, x_m, 2); atol=1e-15)
end

@testset "a sweep without `d` does not fail, it stops checking" begin
    # applicable_relations keeps only relations whose every variable is present,
    # so a bare tuple missing d silently drops the ones that need it. That is
    # the hazard exponents_consistent(::InfiniteRandomness) exists to remove.
    s = InfiniteRandomness(1 // 2, 2 // 1, 1 // 4, 1)
    e = critical_exponents(s)
    with_d = Set(
        nameof(typeof(r)) for r in applicable_relations((; e..., d=s.d); domain=:scaling)
    )
    without_d = Set(nameof(typeof(r)) for r in applicable_relations(e; domain=:scaling))
    @test :ActivatedMomentGrowth in with_d
    @test !(:ActivatedMomentGrowth in without_d)
    @test without_d ⊊ with_d
    @test check_all(e; domain=:scaling)          # green, having checked less
    # the struct door reports on the full set, keyed by relation name
    @test Set(keys(exponent_residuals(s))) ==
        Set(Symbol(lowercase(String(n))) for n in with_d)
end

@testset "inverse map: (ν, ν_typ) recovers ψ at dimension d" begin
    for (ψ, ν, x_m, dd) in _IRFP_SETS
        e = critical_exponents(InfiniteRandomness(ψ, ν, x_m, dd))
        s = infinite_randomness(; ν=e.ν, ν_typ=e.ν_typ, x_m=e.x_m, d=dd)
        @test s.ψ == ψ
        @test critical_exponents(s) == e
    end
end

@testset "the constructor refuses what the name rules out" begin
    @test_throws "not an infinite-randomness fixed point" InfiniteRandomness(
        0, 2, 1 // 4, 1
    )
    @test_throws "ScalingDimensions" InfiniteRandomness(0, 2, 1 // 4, 1)
    @test_throws "not an infinite-randomness fixed point" InfiniteRandomness(
        -1 // 2, 2, 1 // 4, 1
    )
    @test_throws "ν = " InfiniteRandomness(1 // 2, 0, 1 // 4, 1)
    @test_throws "d = " InfiniteRandomness(1 // 2, 2, 1 // 4, 0)
    # ψ = 1 is allowed: ν_typ collapses to 0, which is a statement, not an error
    @test critical_exponents(InfiniteRandomness(1 // 1, 2 // 1, 1 // 4, 1)).ν_typ == 0 // 1
end

@testset "a fixed point that carries `d` needs no second statement of it" begin
    # Both structs answer the registry gate from their own field, which is the
    # only door that cannot be given the wrong dimension.
    for (yt, yh, dd) in _RG_SETS
        s = ScalingDimensions(yt, yh, dd)
        @test exponents_consistent(s)
        @test exponents_consistent(s) == exponents_consistent(critical_exponents(s); d=dd)
        @test exponent_residuals(s) == exponent_residuals(critical_exponents(s); d=dd)
        @test all(iszero, values(exponent_residuals(s)))
    end
    # A wrong dimension is a thing the struct door has no place to accept: the
    # same eigenvalues at another d are a different fixed point, and fail.
    @test !exponents_consistent(
        (critical_exponents(ScalingDimensions(1 // 1, 15 // 8, 2))); d=3
    )
    # and the same for the activated side
    s = InfiniteRandomness(1 // 2, 2 // 1, 1 // 4, 1)
    @test exponent_residuals(s) == exponent_residuals(critical_exponents(s); d=1)
end
