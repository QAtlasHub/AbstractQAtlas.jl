# Bag adoption for the dynamical / spectral verify-engine (issue #77, D1).
#
# The point of externalizing the relation web is that the dynamical /
# spectral consistency network — A(ω), G^R, S(q,ω), χ''(q,ω), their FDT /
# detailed-balance / Dyson ties — is exactly what one cannot hold in one's
# head.  So the turnkey test is: build ONE self-consistent measurement of a
# single damped mode, drop the whole NamedTuple into `relation_report` /
# `check_all`, and confirm every applicable identity fires with correct
# variable-name matching — no per-relation hand-wiring, no pre-projection.
#
# Each pointwise identity is checked against an INDEPENDENT construction
# (the repo's testing contract): A as a Lorentzian vs G^R as a complex
# pole; Dyson's full G assembled from a bare G₀ and Σ; detailed balance
# (an exponential) against the FDT structure factor.

using AbstractQAtlas
using AbstractQAtlas: relation_report, applicable_relations, check_all, check

# names of the relations that fired / are applicable, as a Set of Symbols
_fired(bag; kw...) = Set(nameof(typeof(r.relation)) for r in relation_report(bag; kw...))
_applic(bag; kw...) = Set(nameof(typeof(r)) for r in applicable_relations(bag; kw...))

@testset "one dynamical bag fires the whole pointwise spectral web" begin
    # ---- a self-consistent single damped mode at (q, ω, β) ----
    ω0, η, β = 1.3, 0.08, 1.7
    ω = 0.6                                   # evaluation frequency (single band; q suppressed)

    # retarded Green's function (complex pole) and its spectral weight built
    # the OTHER way, as an explicit Lorentzian — SpectralFromGreens ties them
    GR = 1 / (ω - ω0 + im * η)
    A = (1 / π) * η / ((ω - ω0)^2 + η^2)      # = −Im G^R/π, independent form

    # Dyson: assemble Σ from a bare propagator so the equation closes exactly
    G0 = 1 / (ω - 0.0 + im * η)
    Σ = inv(G0) - inv(GR)

    # antisymmetrized Lorentzian χ''(ω), odd in ω ⇒ FDT gives S(q,ω) and
    # detailed balance S(q,−ω)=e^{−βω}S(q,ω) then holds as a consequence
    _lor(x) = η / (x^2 + η^2) / π
    _χpp(w) = 0.5 * (_lor(w - ω0) - _lor(w + ω0))
    _S(w) = _χpp(w) / (π * (1 - exp(-β * w)))  # DynamicalFDT convention
    χpp, S_plus, S_minus = _χpp(ω), _S(ω), _S(-ω)

    # ---- the measurement bag, in a physicist's natural names ----
    # `q, ω, β, A, G, G0, Σ` feed the spectral representation and Dyson (one
    # complex `G` feeds BOTH — the naming reconciliation that makes this
    # turnkey); `S, χpp` feed the FDT; `S_plus, S_minus` the ±ω detailed
    # balance.  No pre-projection, no per-relation wiring.
    bag = (;
        q=0.0,
        ω=ω,
        β=β,
        A=A,
        G=GR,
        G0=G0,
        Σ=Σ,
        S=S_plus,
        χpp=χpp,
        S_plus=S_plus,
        S_minus=S_minus,
    )

    atol = 1e-10
    fired = _fired(bag; atol=atol, domain=:spectral)

    # every pointwise identity of the dynamical web fires from the one bag …
    @test :Dyson in fired
    @test :SpectralFromGreens in fired
    @test :DynamicalFDT in fired
    @test :DetailedBalance in fired

    # … and every one that fired passes (measurement is self-consistent)
    @test check_all(bag; atol=atol, domain=:spectral)

    # the supplied-integral relations are correctly NOT applicable to a
    # pointwise bag — they need a frequency integral (that is D2 / #19)
    @test :SpectralSumRule ∉ fired
    @test :KramersKronigReal ∉ fired
    @test :StaticFromDynamicalStructureFactor ∉ fired

    # ---- negative control: corrupt one measurement, the web catches it ----
    bad = merge(bag, (; A=2A))                 # wrong spectral weight
    report = relation_report(bad; atol=atol, domain=:spectral)
    failed = Set(nameof(typeof(r.relation)) for r in report if !r.pass)
    @test :SpectralFromGreens in failed        # A no longer matches −Im G^R/π
    @test :DynamicalFDT ∉ failed                # an unrelated identity still holds
    @test !check_all(bad; atol=atol, domain=:spectral)
end

@testset "supplied-integral relations fire once the integral field is present" begin
    # A single retarded pole: its real/imag parts form a Kramers–Kronig pair
    # and its spectral weight integrates to 1 — closed forms we can supply
    # directly here.  D2 (#19) is what COMPUTES these integrals from a
    # spectrum (grid or pole–residue); this test only confirms the relations
    # fire and match variable names once the integral value is in the bag.
    ω0, η, ω = 0.7, 0.1, 0.35
    GR = 1 / (ω - ω0 + im * η)

    # the pointwise bag alone does not reach these relations …
    pointwise = (; ω=ω, A=-imag(GR) / π, G=GR)
    @test :SpectralSumRule ∉ _applic(pointwise; domain=:spectral)
    @test :KramersKronigReal ∉ _applic(pointwise; domain=:spectral)

    # … but adding the (here analytic) integrals makes them applicable + pass:
    #   ∫A dω = 1 (Lorentzian normalization)
    #   Re G^R(ω) = (1/π) P∫ Im G^R(ω')/(ω'−ω) dω'   ⇒  pv_imag = π Re G^R
    #   Im G^R(ω) = −(1/π) P∫ Re G^R(ω')/(ω'−ω) dω'  ⇒  pv_real = −π Im G^R
    enriched = merge(
        pointwise,
        (;
            spectral_integral=1.0,
            Reχ=real(GR),
            pv_imag=π * real(GR),
            Imχ=imag(GR),
            pv_real=-π * imag(GR),
        ),
    )
    applic = _applic(enriched; domain=:spectral)
    @test :SpectralSumRule in applic
    @test :KramersKronigReal in applic
    @test :KramersKronigImag in applic

    @test check(SpectralSumRule(); spectral_integral=1.0, atol=0)
    @test check_all(enriched; atol=1e-10, domain=:spectral)
end
