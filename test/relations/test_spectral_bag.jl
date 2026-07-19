# Bag adoption for the dynamical / spectral verify-engine (issue #77, item 1).
#
# The point of externalizing the relation web is that the dynamical / spectral
# consistency network — A(ω), G^R, S(q,ω), χ''(q,ω), their FDT / detailed-balance /
# Dyson ties — is exactly what one cannot hold in one's head.  So the turnkey test is:
# build ONE self-consistent measurement of a single damped mode, drop the whole
# NamedTuple into `relation_report` / `check_all`, and confirm every applicable
# identity fires with correct variable-name matching — no per-relation hand-wiring, no
# pre-projection (one complex `G` feeds BOTH Dyson and the spectral representation).
#
# Each pointwise identity is checked against an INDEPENDENT construction (the repo's
# testing contract): A as a Lorentzian vs G^R as a complex pole; Dyson's Σ assembled
# from a bare G₀; detailed balance (an exponential) as a consequence of the FDT.

using AbstractQAtlas
using AbstractQAtlas: relation_report, applicable_relations, check_all, check

# names of the relations that fired / are applicable, as a Set of Symbols
_fired(bag; kw...) = Set(nameof(typeof(r.relation)) for r in relation_report(bag; kw...))
_applic(bag; kw...) = Set(nameof(typeof(r)) for r in applicable_relations(bag; kw...))

@testset "one dynamical bag fires the whole pointwise spectral web" begin
    # ---- a self-consistent single damped mode at (q, ω, β) ----
    ω0, η, β = 1.3, 0.08, 1.7
    ω = 0.6                                   # evaluation frequency (single band; q suppressed)

    # retarded Green's function (complex pole) and its spectral weight built the OTHER
    # way, as an explicit Lorentzian — SpectralFromGreens ties them
    GR = 1 / (ω - ω0 + im * η)
    A = (1 / π) * η / ((ω - ω0)^2 + η^2)      # = −Im G^R/π, independent form

    # Dyson: assemble Σ from a bare propagator so the equation closes exactly
    G0 = 1 / (ω - 0.0 + im * η)
    Σ = inv(G0) - inv(GR)

    # antisymmetrized Lorentzian χ''(ω), odd in ω ⇒ FDT gives S(q,ω) and detailed
    # balance S(q,−ω) = e^{−βω} S(q,ω) then holds as a consequence
    _lor(x) = η / (x^2 + η^2) / π
    _χpp(w) = 0.5 * (_lor(w - ω0) - _lor(w + ω0))
    _S(w) = _χpp(w) / (π * (1 - exp(-β * w)))  # DynamicalFDT convention
    χpp, S_plus, S_minus = _χpp(ω), _S(ω), _S(-ω)

    # ---- the measurement bag, in a physicist's natural names ----
    # one complex `G` feeds BOTH Dyson and the spectral representation (the naming
    # reconciliation that makes this turnkey); `S, χpp` feed the FDT; `S_plus, S_minus`
    # the ±ω detailed balance.  No pre-projection, no per-relation wiring.
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

    # … and every one that fired passes (the measurement is self-consistent)
    @test check_all(bag; atol=atol, domain=:spectral)

    # the supplied-integral relations are correctly NOT applicable to a pointwise bag —
    # they need a frequency integral (that auto-evaluation is item 2 / #19)
    @test :SpectralSumRule ∉ fired
    @test :KramersKronigReal ∉ fired
    @test :StaticFromDynamicalStructureFactor ∉ fired
    @test :FSumRule ∉ fired

    # ---- negative control: corrupt one measurement, the web catches it ----
    bad = merge(bag, (; A=2A))                 # wrong spectral weight
    report = relation_report(bad; atol=atol, domain=:spectral)
    failed = Set(nameof(typeof(r.relation)) for r in report if !r.pass)
    @test :SpectralFromGreens in failed        # A no longer matches −Im G^R/π
    @test :DynamicalFDT ∉ failed                # an unrelated identity still holds
    @test !check_all(bad; atol=atol, domain=:spectral)
end

@testset "supplied-integral relations fire once the caller adds the integrals" begin
    # the same web, now WITH the caller-computed frequency integrals: the sum-rule /
    # static-SF / f-sum relations that were (correctly) inapplicable pointwise now fire
    # too — the supplied-integral convention makes the whole web gateable WITHOUT the
    # package computing any transform (that auto-evaluation from a grid is #19).
    q, m = 2.0, 1.0
    bag = (;
        spectral_integral=1.0,                 # ∫A dω = 1
        Sq=0.5,
        sqw_integral=0.5,                       # S(q) = ∫S(q,ω)dω
        first_moment=q^2 / (2 * m),             # ∫ω S(q,ω)dω = q²/2m (N=1)
        q=q,
        m=m,
    )
    fired = _fired(bag; atol=1e-12, domain=:spectral)
    @test :SpectralSumRule in fired
    @test :StaticFromDynamicalStructureFactor in fired
    @test :FSumRule in fired
    @test check_all(bag; atol=1e-12, domain=:spectral)
end

@testset "multi-time χ⁽²⁾ bag fires the reality web (#77 item 3)" begin
    # a real-in-time 2nd-order response χ⁽²⁾(ω1,ω2)=∏ 1/(ω0²−ωᵢ²−iγωᵢ): reality holds, so
    # measuring it at ω⃗ and −ω⃗ and dropping both in the bag fires the Re-even/Im-odd web
    ω0, γ = 1.4, 0.1
    χ2(w1, w2) = 1 / ((ω0^2 - w1^2 - im * γ * w1) * (ω0^2 - w2^2 - im * γ * w2))
    cp, cm = χ2(0.4, 0.9), χ2(-0.4, -0.9)      # χ⁽²⁾(ω⃗), χ⁽²⁾(−ω⃗)
    bag = (; Re_plus=real(cp), Re_minus=real(cm), Im_plus=imag(cp), Im_minus=imag(cm))
    fired = _fired(bag; atol=1e-10, domain=:spectral)
    @test :ResponseRealityReal in fired && :ResponseRealityImag in fired
    @test check_all(bag; atol=1e-10, domain=:spectral)
end
