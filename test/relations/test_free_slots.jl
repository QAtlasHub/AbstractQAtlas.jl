# A slot left free absorbs the data, and the relation then cannot fail.
#
# Each case below was MEASURED passing before the fix, on inputs that violate the
# law by orders of magnitude. The pattern is the registry's own: a value the
# physics fixes (a universal constant, a two-valued label, an equilibrium
# distribution) written as a caller-supplied slot with nothing pinning it.

using AbstractQAtlas
using Test: @test, @test_throws, @testset

const L_SOMMERFELD = π^2 / 3

why(f) =
    try
        f()
        ""
    catch e
        sprint(showerror, e)
    end

@testset "the Lorenz number is a constant, so omitting it tests the law" begin
    σ, T = 2.0, 300.0
    obeys = L_SOMMERFELD * σ * T
    # Omitted: the default IS the law, so a material off by 1000 fails.
    @test check(WiedemannFranz(); κ=obeys, σ=σ, T=T)
    @test !check(WiedemannFranz(); κ=1000 * obeys, σ=σ, T=T)
    @test !check(WiedemannFranz(); κ=(-obeys), σ=σ, T=T)   # negative κ is not a metal
    # Supplied: still accepted, because a non-Fermi liquid has its own Lorenz
    # number. The docstring says so, and says what a self-derived L0 is worth.
    @test check(WiedemannFranz(); κ=1000 * obeys, σ=σ, T=T, L0=1000 * L_SOMMERFELD)
    # Solving FOR it is the honest way to get the ratio out of a measurement.
    @test solve(WiedemannFranz(), Val(:L0); κ=obeys, σ=σ, T=T) ≈ L_SOMMERFELD

    xy = L_SOMMERFELD * T * 0.5
    @test check(RighiLeduc(); κxy=xy, T=T, σxy=0.5)
    @test !check(RighiLeduc(); κxy=1e6 * xy, T=T, σxy=0.5)
end

@testset "the statistics sign takes two values, so a third is refused" begin
    β, ω, Ggtr = 1.0, 2.0, 1.0
    kms(ζ) =
        check(KMSGreaterLesser(); Gles=ζ * exp(-β * ω) * Ggtr, Ggtr=Ggtr, ζ=ζ, ω=ω, β=β)
    @test kms(+1)
    @test kms(-1)
    # Free, ζ absorbs any ratio: `ζ = G^</(e^{-βω} G^>)` zeroes the residual for
    # every pair, so the KMS condition could never fail. 7.3 passed before this.
    @test occursin("exchange-statistics sign", why(() -> kms(7.3)))
    @test occursin("pass without being tested", why(() -> kms(0)))
    # And the condition itself still bites: a pair off the KMS ratio fails at a
    # legal ζ, which is what says the guard did not simply replace one tautology
    # with another.
    @test !check(
        KMSGreaterLesser(); Gles=5.0 * exp(-β * ω) * Ggtr, Ggtr=Ggtr, ζ=1, ω=ω, β=β
    )
end

@testset "the FDT needs the half that pins h" begin
    β, ω = 1.0, 2.0
    GR, GA = 1.0 + 2.0im, 1.0 - 2.0im
    for stat in (Bosonic(), Fermionic())
        h_eq = keldysh_distribution(stat, ω; β=β)
        # `KeldyshFDT` reads h as given, so it holds for ANY state. That is not a
        # defect to fix in it: `G^K = h(G^R-G^A)` is the definition of h.
        for h in (h_eq, 0.017, -h_eq)
            @test check(KeldyshFDT(); GK=h * (GR - GA), h=h, GR=GR, GA=GA)
        end
        # The claim "this system is thermal" lives in the other half, and it fails.
        eq(h) =
            check(KeldyshDistributionEquilibrium(); h=h, ω=ω, β=β, stat=stat, atol=1e-10)
        @test eq(h_eq)
        @test !eq(0.017)
        @test !eq(-h_eq)
        # The self-energy mirror shares the same h, so one relation covers both.
        ΣR, ΣA = 0.3 + 0.7im, 0.3 - 0.7im
        @test check(
            SelfEnergyKeldyshFDT(); SigmaK=h_eq * (ΣR - ΣA), h=h_eq, SigmaR=ΣR, SigmaA=ΣA
        )
    end
    # Bosonic and fermionic h differ, so declaring the wrong statistics is caught.
    @test !check(
        KeldyshDistributionEquilibrium();
        h=keldysh_distribution(Bosonic(), ω; β=β),
        ω=ω,
        β=β,
        stat=Fermionic(),
        atol=1e-10,
    )
end
