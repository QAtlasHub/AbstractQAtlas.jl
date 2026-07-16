# The Keldysh RAK structure and its equilibrium fluctuation–dissipation
# constraint.  The RAK relations are algebraic identities (exact); the headline
# check is that the fluctuation–dissipation theorem is a CONSEQUENCE of the KMS
# detailed-balance relation plus those identities — verified against the
# INDEPENDENT distribution function and the existing spectral definition, not
# against the Keldysh relations' own output.

using AbstractQAtlas
using AbstractQAtlas:
    residual,
    check,
    solve,
    KeldyshComponent,
    KeldyshCausality,
    AdvancedRetardedConjugate,
    KeldyshFDT,
    KMSGreaterLesser,
    SpectralFromKeldysh,
    SpectralFromGreens,
    keldysh_distribution,
    occupation,
    Fermionic,
    Bosonic,
    spectral_origin,
    origin_relation,
    operation_scope,
    all_relations

@testset "RAK identities are exact (equilibrium or not)" begin
    # arbitrary greater/lesser correlators — the component/causality relations
    # are pure definitions, exact in rationals
    Ggtr, Gles = 5 // 2, -3 // 4
    GK = Ggtr + Gles
    GRmGA = Ggtr - Gles
    @test residual(KeldyshComponent(); GK=GK, Ggtr=Ggtr, Gles=Gles) == 0 // 1
    # KeldyshCausality with any GR, GA whose difference matches
    GR, GA = 1 // 1, 1 - GRmGA          # GR − GA = GRmGA
    @test residual(KeldyshCausality(); GR=GR, GA=GA, Ggtr=Ggtr, Gles=Gles) == 0 // 1
    # advanced = conjugate of retarded (complex)
    z = 0.3 - 0.8im
    @test residual(AdvancedRetardedConjugate(); GA=conj(z), GR=z) == 0
    @test !check(AdvancedRetardedConjugate(); GA=z, GR=z)   # z ≠ conj(z) for complex z
end

@testset "FDT is a CONSEQUENCE of KMS + the RAK identities" begin
    # Give ONLY the greater correlator and the KMS/detailed-balance law
    # G^< = ζ e^{−βω} G^>.  Then G^K/(G^R−G^A) must equal the independent
    # distribution function, and KeldyshFDT must hold — for both statistics.
    for (stat, ζ) in ((Bosonic(), 1), (Fermionic(), -1)), β in (0.5, 1.3), ω in (0.4, 2.1)
        Ggtr = 1.7 - 0.9im
        Gles = ζ * exp(-β * ω) * Ggtr
        # KMS relation is satisfied by construction
        @test residual(KMSGreaterLesser(); Gles=Gles, Ggtr=Ggtr, ζ=ζ, ω=ω, β=β) ≈ 0 atol =
            1e-12
        GRmGA = Ggtr - Gles                 # = G^R − G^A  (KeldyshCausality)
        GK = Ggtr + Gles                    # = G^K        (KeldyshComponent)
        h = keldysh_distribution(stat, ω; β=β)
        # the ratio G^K/(G^R−G^A) IS the distribution function (the theorem)
        @test GK / GRmGA ≈ h
        # …hence KeldyshFDT holds with that h (GA absorbed into GRmGA via GA=0)
        @test residual(KeldyshFDT(); GK=GK, h=h, GR=GRmGA, GA=0) ≈ 0 atol = 1e-12
        # accept T in place of β at the verb layer
        @test keldysh_distribution(stat, ω; T=1 / β) ≈ h
    end
end

@testset "keldysh_distribution = 1 ∓ 2n and its analytic structure" begin
    for β in (0.7, 1.9), ω in (0.3, 1.4)
        hF = keldysh_distribution(Fermionic(), ω; β=β)
        hB = keldysh_distribution(Bosonic(), ω; β=β)
        # tie to the Fermi–Dirac / Bose–Einstein occupation (μ = 0)
        @test hF ≈ 1 - 2 * occupation(Fermionic(), ω; β=β)
        @test hB ≈ 1 + 2 * occupation(Bosonic(), ω; β=β)
        # h is odd in ω
        @test keldysh_distribution(Fermionic(), -ω; β=β) ≈ -hF
        @test keldysh_distribution(Bosonic(), -ω; β=β) ≈ -hB
        # fermionic h is bounded in [-1, 1]; bosonic |h| ≥ 1
        @test abs(hF) ≤ 1
        @test abs(hB) ≥ 1
    end
end

@testset "SpectralFromKeldysh bridges to A = −Im Gᴿ/π" begin
    # a single-pole retarded propagator; A = −Im Gᴿ/π by SpectralFromGreens
    for (ε, η) in ((0.4, 0.05), (-1.2, 0.1)), ω in (-0.3, 0.9)
        GR = 1 / (ω - ε + im * η)
        GA = conj(GR)                       # AdvancedRetardedConjugate
        A = -imag(GR) / π
        # Keldysh spectral bridge reproduces the same A …
        @test residual(SpectralFromKeldysh(); A=A, GR=GR, GA=GA) ≈ 0 atol = 1e-12
        # … and agrees with the existing Matsubara-side definition
        @test residual(SpectralFromGreens(); A=A, G=GR) ≈ 0 atol = 1e-12
    end
end

@testset "generic solve inverts every Keldysh variable it is affine in" begin
    Ggtr, Gles = 2.0 + 0.1im, -0.5 + 0.2im
    GK, GRmGA = Ggtr + Gles, Ggtr - Gles
    # KeldyshComponent solved for each correlator
    @test solve(KeldyshComponent(), Val(:GK); Ggtr=Ggtr, Gles=Gles) ≈ GK
    @test solve(KeldyshComponent(), Val(:Gles); GK=GK, Ggtr=Ggtr) ≈ Gles
    # KeldyshFDT solved for the Keldysh component and the distribution function
    h = 1.3
    @test solve(KeldyshFDT(), Val(:GK); h=h, GR=GRmGA, GA=0) ≈ h * GRmGA
    @test solve(KeldyshFDT(), Val(:h); GK=h * GRmGA, GR=GRmGA, GA=0) ≈ h
end

@testset "G^A joins the spectral graph via the adjoint edge" begin
    o = spectral_origin(AdvancedGreensFunction)
    @test o !== nothing
    @test o.from === RetardedGreensFunction
    @test o.via === :adjoint
    # the adjoint edge is pointwise (definitional), tied to its relation
    @test origin_relation(:adjoint) isa AdvancedRetardedConjugate
    @test operation_scope(:adjoint) === :definitional
end

@testset "Keldysh relations register under :keldysh" begin
    @test length(all_relations(; domain=:keldysh)) == 6
end
