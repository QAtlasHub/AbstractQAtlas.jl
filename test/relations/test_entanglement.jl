# Entanglement-entropy relations vs INDEPENDENT constructions: purity from an
# explicit density matrix, the cut-counting coefficient from the EXACT
# free-fermion critical Ising chain, and Page's formula against exact small cases
# and a Haar-random-state average.  The chain is quadratic in Majoranas, so END
# and CENTRED blocks come from the SAME exact ground state — the cut count is the
# only difference between them.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve
using LinearAlgebra, Random

using ExperimentalAPI: ExperimentalAPI

struct _UnknownBC <: AbstractQAtlas.BoundaryCondition end

@testset "Rényi-2 from purity: S_2 = −ln Tr ρ²" begin
    # pure state: purity 1, S_2 = 0
    @test solve(RenyiTwoPurity(), Val(:S2); purity=1.0) == 0.0
    # maximally mixed on d levels: purity = 1/d, S_2 = ln d
    for d in (2, 4, 8)
        ρ = Matrix{Float64}(I, d, d) ./ d
        purity = tr(ρ^2)                      # = 1/d, built independently
        @test check(RenyiTwoPurity(); S2=log(d), purity=purity, atol=1e-13)
    end
    # a generic ρ: S_2 vs −ln Tr ρ² from the explicit matrix
    ρ = [0.6 0.1; 0.1 0.4]
    @test check(RenyiTwoPurity(); S2=(-log(tr(ρ^2))), purity=tr(ρ^2), atol=1e-13)
end

@testset "CFT entanglement slope counts cuts: exact critical Ising chain" begin
    # CLEAN critical TFIM, open ends, c = 1/2: H = -Σ Z_j Z_{j+1} - Σ X_j is
    # (i/4) Σ A_{mn} a_m a_n with A[2j-1,2j] = -2h_j, A[2j,2j+1] = -2J_j.  The
    # ground-state covariance is A's orthogonal polar factor; a region's entropy
    # is Peschel's function of its restricted spectrum.
    c = 1 / 2

    function covariance(N)
        A = zeros(2N, 2N)
        for j in 1:N
            A[2j - 1, 2j] = -2.0        # -2 h_j, h_j = 1
        end
        for j in 1:(N - 1)
            A[2j, 2j + 1] = -2.0        # -2 J_j, J_j = 1
        end
        A = A - transpose(A)            # antisymmetry as structure, not as bookkeeping
        F = svd(A)
        Γ = F.U * F.Vt
        return (Γ .- transpose(Γ)) ./ 2   # antisymmetry drifts at a zero singular value
    end

    function entropy(Γ, sites)
        idx = vcat(([2j - 1, 2j] for j in sites)...)
        ν = eigvals(Hermitian(im .* Γ[idx, idx]))
        return free_fermion_entanglement_entropy([
            (1 + real(x)) / 2 for x in ν if real(x) > 0
        ])
    end

    # OLS slope of S against ln ℓ; the additive constant drops out of the weights
    function log_slope(ℓs, S)
        x = log.(float.(ℓs))
        x .-= sum(x) / length(x)
        return sum(x .* S) / sum(abs2, x)
    end

    # ONE cut: the block sits at an open end, so only its inner edge is a cut.
    # TWO cuts: the block sits in the bulk, both edges away from the boundary.
    end_block(N, ℓ) = 1:ℓ
    centred_block(N, ℓ) = (div(N - ℓ, 2) + 1):(div(N - ℓ, 2) + ℓ)

    measured = map((48, 128)) do N
        Γ = covariance(N)
        ℓs = 3:div(N, 4)          # short of ℓ ≈ N, where a centred block fills the chain
        a1 = log_slope(ℓs, [entropy(Γ, end_block(N, ℓ)) for ℓ in ℓs])
        a2 = log_slope(ℓs, [entropy(Γ, centred_block(N, ℓ)) for ℓ in ℓs])
        (; N, a1, a2)
    end

    for m in measured
        # Budgets are the MEASURED deviation rounded up; the fixture is
        # deterministic to ~1e-14, so they are not noise budgets:
        #
        #   ncuts=1:  0.000166 (N=48), 0.000179 (N=128)  →  atol 0.0005  (2.8×)
        #   ncuts=2:  0.011596 (N=48), 0.007269 (N=128)  →  atol 0.015   (1.3×)
        #
        # Two cuts carries the larger correction — both edges sit a finite
        # distance from the chain's ends — and it shrinks with N, as it must.
        @test check(CFTEntanglementSlope(); dS_dlogℓ=m.a1, c=c, ncuts=1, atol=0.0005)
        @test check(CFTEntanglementSlope(); dS_dlogℓ=m.a2, c=c, ncuts=2, atol=0.015)

        # ...and it can DISAGREE: each geometry is nearer its own cut count than
        # the other's, which a relation ignoring `ncuts` could not satisfy.
        @test abs(m.a1 - 1 * c / 6) < abs(m.a1 - 2 * c / 6)
        @test abs(m.a2 - 2 * c / 6) < abs(m.a2 - 1 * c / 6)

        # `c` cancels, so the ratio is `ncuts` alone.  Measured |ratio−2| = 0.143
        # (N=48), 0.091 (N=128), same ~1.3× round-up; one cut would put it at 1.
        @test m.a2 / m.a1 ≈ 2 atol = 0.18

        # The same measured slope yields a DIFFERENT central charge under a
        # different cut count.
        @test solve(CFTEntanglementSlope(), Val(:c); dS_dlogℓ=m.a2, ncuts=2) ≈ c rtol = 0.08
        @test solve(CFTEntanglementSlope(), Val(:c); dS_dlogℓ=m.a2, ncuts=1) ≈ 2c rtol =
            0.08
    end

    # An axis that must not matter: the coefficient belongs to the fixed point,
    # so a larger chain moves the ratio TOWARDS 2, never away.
    @test abs(measured[2].a2 / measured[2].a1 - 2) <
        abs(measured[1].a2 / measured[1].a1 - 2)

    # Pure arithmetic on the declared form, with no fixture in the way.
    @test check(CFTEntanglementSlope(); dS_dlogℓ=1 / 3, c=1.0, ncuts=2, atol=1e-14)
    @test check(CFTEntanglementSlope(); dS_dlogℓ=1 / 6, c=1.0, ncuts=1, atol=1e-14)

    # Data that does not say which geometry the slope was measured on yields NO
    # central-charge row, rather than one that assumed a geometry for the caller.
    @test isempty(relation_report((; dS_dlogℓ=1 / 3, c=1.0)))
    discovered = relation_report((; dS_dlogℓ=1 / 3, c=1.0, ncuts=2))
    @test only(discovered).relation isa CFTEntanglementSlope
    @test only(discovered).pass
end

@testset "InfiniteRandomnessEntanglementSlope reproduces Refael-Moore" begin
    # Refael & Moore 2004: c̃ = (ln 2)/2 for the random transverse-field Ising
    # chain (Eq. 23) and ln 2 for the random singlet phase of the Heisenberg and
    # XX chains (Eq. 19).
    c̃_ising, c̃_singlet = log(2) / 2, log(2)

    # Their published slopes, for a segment of a chain, which is two cuts.
    @test check(
        InfiniteRandomnessEntanglementSlope();
        dS_dlogℓ=log(2) / 6,
        c̃=c̃_ising,
        ncuts=2,
        atol=1e-14,
    )
    @test check(
        InfiniteRandomnessEntanglementSlope();
        dS_dlogℓ=log(2) / 3,
        c̃=c̃_singlet,
        ncuts=2,
        atol=1e-14,
    )

    # Eqs. (19) and (23) are written in bits, `S = -Tr ρ log₂ ρ`, where the same
    # two-cut slopes read 1/3 and 1/6 exactly.  Dividing by ln 2 must land there,
    # or the relation is off by the base.
    @test (log(2) / 3) / log(2) ≈ 1 / 3 atol = 1e-14
    @test (log(2) / 6) / log(2) ≈ 1 / 6 atol = 1e-14

    # Getting the cut count wrong is not a small error: read at one cut, the
    # random Ising slope returns the random singlet value, so the two classes
    # trade places.  This is the failure the `ncuts` axis exists to prevent.
    @test solve(
        InfiniteRandomnessEntanglementSlope(), Val(:c̃); dS_dlogℓ=log(2) / 6, ncuts=2
    ) ≈ c̃_ising atol = 1e-14
    @test solve(
        InfiniteRandomnessEntanglementSlope(), Val(:c̃); dS_dlogℓ=log(2) / 6, ncuts=1
    ) ≈ c̃_singlet atol = 1e-14

    # A block at an open end is one cut: half the slope, ln 2/12 = 0.0578, not
    # the ln 2/6 = 0.1155 of a segment.  The two must not be interchangeable.
    @test check(
        InfiniteRandomnessEntanglementSlope();
        dS_dlogℓ=log(2) / 12,
        c̃=c̃_ising,
        ncuts=1,
        atol=1e-14,
    )
    @test !check(
        InfiniteRandomnessEntanglementSlope();
        dS_dlogℓ=log(2) / 12,
        c̃=c̃_ising,
        ncuts=2,
        atol=1e-3,
    )

    # As in the CFT case, a slope that does not say its geometry yields no row.
    @test isempty(relation_report((; dS_dlogℓ=log(2) / 6, c̃=c̃_ising)))
    found = relation_report((; dS_dlogℓ=log(2) / 6, c̃=c̃_ising, ncuts=2))
    @test only(found).relation isa InfiniteRandomnessEntanglementSlope
    @test only(found).pass
end

@testset "entanglement_cuts derives b from the Region and the boundary condition" begin
    # A contiguous block has two cuts on a ring, and one only when it reaches an
    # end of an open chain. The same sites give different answers, which is the
    # whole reason the count cannot be read off the boundary condition.
    @test entanglement_cuts(PBC(8), Region(2, 3, 4)) == 2
    @test entanglement_cuts(OBC(8), Region(2, 3, 4)) == 2
    @test entanglement_cuts(OBC(8), Region(1, 2, 3)) == 1
    @test entanglement_cuts(OBC(8), Region(6, 7, 8)) == 1

    # Filling the system leaves nothing to cut, on either boundary condition.
    @test entanglement_cuts(PBC(4), Region(1, 2, 3, 4)) == 0
    @test entanglement_cuts(OBC(4), Region(1, 2, 3, 4)) == 0

    # The ring's wrap-around is a real adjacency: a block straddling it stays at
    # two cuts, where an open chain would count the same sites as two blocks.
    @test entanglement_cuts(PBC(8), Region(8, 1, 2)) == 2
    @test entanglement_cuts(OBC(8), Region(8, 1, 2)) == 2
    @test entanglement_cuts(PBC(8), Region(2, 3, 6, 7)) == 4

    # An infinite chain needs no N: adjacency around the region is enough.
    @test entanglement_cuts(Infinite(), Region(5, 6, 7)) == 2
    @test entanglement_cuts(Infinite(), Region(5, 6, 9)) == 4
    @test entanglement_cuts(PBC(8), Region()) == 0

    # A set with no adjacency is refused rather than guessed at, and so are sites
    # the declared chain does not contain, at BOTH ends: a bare `ErrorException`
    # would not say the intended guard is the one that fired, and reordering them
    # would swap the diagnoses without failing.
    @test_throws "adjacency needs integer sites" entanglement_cuts(OBC(8), Region("a", "b"))
    @test_throws "fall outside the chain" entanglement_cuts(OBC(4), Region(3, 4, 5))
    @test_throws "fall outside the chain" entanglement_cuts(OBC(4), Region(0, 1, 2))
    @test_throws "declares no chain length" entanglement_cuts(OBC(), Region(1, 2))
    @test_throws "declares no chain length" entanglement_cuts(PBC(), Region(1, 2))
end

@testset "finite-size entropy: ring, open chain, and what separates them" begin
    c, c₁, L = 1 / 2, 0.4785, 2048.0          # Ising; c₁ non-universal, cancels below

    ring(ℓ, Lc=L) = (c / 3) * log((Lc / π) * sin(π * ℓ / Lc)) + c₁
    open_(ℓ, Lc=L) = (c / 6) * log((2Lc / π) * sin(π * ℓ / Lc)) + 0.0 + c₁ / 2

    @test check(CFTEntanglementPBC(); S=ring(512), c=c, L=L, ℓ=512, c₁=c₁, atol=1e-12)
    # Iglói & Lin Table 1 measure g → 1 for the Ising chain, so `ln g` vanishes there.
    @test check(
        CFTEntanglementOBC(); S=open_(512), c=c, L=L, ℓ=512, c₁=c₁, ln_g=0.0, atol=1e-12
    )

    # A ring block far from filling the ring is the infinite chain: the chord
    # tends to ℓ, so the two forms must agree, and increasingly so.
    inf_chain(ℓ) = (c / 3) * log(ℓ) + c₁
    @test abs(ring(8) - inf_chain(8)) < abs(ring(256) - inf_chain(256))
    @test ring(8) ≈ inf_chain(8) atol = 2e-5

    # The open chain is NOT the ring halved. Same (L, ℓ, c, c₁), and the gap is
    # the chord's 2L/π against L/π; assuming a bare factor of two misses it.
    @test !isapprox(open_(512), ring(512) / 2; atol=1e-6)
    @test open_(512) - ring(512) / 2 ≈ (c / 6) * log(2) atol = 1e-12
end

@testset "HalvedChainEntropyDifference is exact on both boundary conditions" begin
    c, c₁, L = 1 / 2, 0.4785, 4096.0
    ring(ℓ, Lc) = (c / 3) * log((Lc / π) * sin(π * ℓ / Lc)) + c₁
    open_(ℓ, Lc) = (c / 6) * log((2Lc / π) * sin(π * ℓ / Lc)) + c₁ / 2

    # Derived from the two forms above, not assumed: halving shifts the chord by
    # exactly two, and c₁ (and `ln g`) cancel, which is why no constant appears.
    ΔS_ring = ring(L / 2, L) - ring(L / 4, L / 2)
    ΔS_open = open_(L / 2, L) - open_(L / 4, L / 2)
    @test check(HalvedChainEntropyDifference(); ΔS=ΔS_ring, c=c, ncuts=2, atol=1e-12)
    @test check(HalvedChainEntropyDifference(); ΔS=ΔS_open, c=c, ncuts=1, atol=1e-12)

    # The `ln 2` is load-bearing, and its absence is not merely a wrong number.
    # The source reads ΔS = c/3 and c/6 because it counts bits; applying that to
    # entropies in nats returns c·ln 2 rather than c.
    @test solve(HalvedChainEntropyDifference(), Val(:c); ΔS=ΔS_ring, ncuts=2) ≈ c atol =
        1e-12
    @test ΔS_ring / (2 / 6) ≈ c * log(2) atol = 1e-12
    @test !isapprox(ΔS_ring / (2 / 6), c; atol=1e-3)

    # For the Ising chain that misreading lands exactly on the random chain's
    # effective central charge, since c̃ = c ln 2 there: a clean chain measured
    # in the wrong base is numerically an infinite-randomness one.
    @test c * log(2) ≈ log(2) / 2 atol = 1e-15

    # And the cut count still discriminates: the same difference read at the
    # wrong count returns twice or half the central charge.
    @test solve(HalvedChainEntropyDifference(), Val(:c); ΔS=ΔS_ring, ncuts=1) ≈ 2c atol =
        1e-12
end

@testset "off-critical saturation counts the same boundary points" begin
    c, ξ = 1 / 2, 40.0
    # Iglói & Lin write the prefactor as `b`, the number of boundary points.
    @test check(
        OffCriticalEntanglementSaturation();
        S=2 * (c / 6) * log(ξ),
        c=c,
        ξ=ξ,
        ncuts=2,
        atol=1e-12,
    )
    @test check(
        OffCriticalEntanglementSaturation();
        S=1 * (c / 6) * log(ξ),
        c=c,
        ξ=ξ,
        ncuts=1,
        atol=1e-12,
    )
    # ξ replaces ℓ: at ξ = ℓ the saturated value meets the critical logarithm.
    @test 2 * (c / 6) * log(ξ) ≈ (c / 3) * log(ξ) atol = 1e-12
end

@testset "the guards refuse what the raw relations cannot see" begin
    # `sin` is periodic, so ℓ outside the chain aliases onto a legitimate answer
    # rather than looking wrong: ℓ = 250 on L = 100 returned exactly the ℓ = 50
    # value before the guard.
    @test_throws "need 0 < ℓ < L" AbstractQAtlas.solve(
        CFTEntanglementPBC(), Val(:S); c=0.5, L=100.0, ℓ=250.0, c₁=0.4785
    )
    @test AbstractQAtlas.solve(
        CFTEntanglementPBC(), Val(:S); c=0.5, L=100.0, ℓ=50.0, c₁=0.4785
    ) isa Real

    # ℓ = L is not caught by a blow-up: `sin(float(π))` is 1.2e-16, so the answer
    # was a large finite number set by rounding, where the truth is 0.
    @test_throws "need 0 < ℓ < L" AbstractQAtlas.solve(
        CFTEntanglementOBC(), Val(:S); c=0.5, L=100.0, ℓ=100.0, c₁=0.4785, ln_g=0.0
    )

    # Solving for ℓ is refused, and by the guard rather than by affinity: `_solve`
    # probes at ℓ = 0 first, so the domain check fires before the parabola test is
    # reached. Pinned to the message, since a bare Exception cannot tell which.
    @test_throws "need 0 < ℓ < L" AbstractQAtlas.solve(
        CFTEntanglementPBC(), Val(:ℓ); S=1.0, c=0.5, L=100.0, c₁=0.4
    )

    # Unsigned labels: `lo - 1` wraps to typemax and empties the range, which
    # reported 0 cuts for a two-cut region.
    @test entanglement_cuts(Infinite(), Region(UInt(0), UInt(1))) ==
        entanglement_cuts(Infinite(), Region(0, 1)) ==
        2
    @test_throws "not a lattice index" entanglement_cuts(OBC(8), Region(true))

    # A boundary condition with no branch must not read as an open chain.
    @test_throws "no adjacency defined" entanglement_cuts(_UnknownBC(), Region(1, 2))

    # The random ring's `f` is only a number to the relation, so the sign guard is
    # all it can offer; the geometry-aware sampling lives in the region sweep.
    @test_throws "L = -1024" AbstractQAtlas.solve(
        InfiniteRandomnessEntanglementPBC(), Val(:S̄); c̃=0.25, L=-1024.0, f=-0.3, c₁′=0.3
    )
end

@testset "only the unsettled half is marked experimental" begin
    _QI = AbstractQAtlas.QuantumInformation

    # Eq. (24) has no independent oracle here and `f` admits no bound, so the two
    # names that carry it say so at runtime.
    @test ExperimentalAPI.isexperimental(_QI, :InfiniteRandomnessEntanglementPBC)

    # The rest must NOT be marked. A flag on everything reports nothing, and these
    # are anchored to published constants with discriminating tests above.
    @test !ExperimentalAPI.isexperimental(_QI, :InfiniteRandomnessEntanglementSlope)
    @test !ExperimentalAPI.isexperimental(_QI, :CFTEntanglementPBC)
    @test !ExperimentalAPI.isexperimental(_QI, :CFTEntanglementOBC)
    @test !ExperimentalAPI.isexperimental(_QI, :entanglement_cuts)
end

@testset "the conformal chord is the one-harmonic case of the random one" begin
    c̃, c₁′, L, ℓ = log(2) / 2, 0.31, 1024.0, 300.0

    # `Σₖ Aₖ(2k-1)π = 1` with a single harmonic forces A₁ = 1/π, and then
    # `L f(ℓ/L)` is the conformal chord exactly.
    A₁ = 1 / π
    @test A₁ * (2 * 1 - 1) * π ≈ 1 atol = 1e-14
    f_conformal(v) = A₁ * sin(π * v)
    @test L * f_conformal(ℓ / L) ≈ (L / π) * sin(π * ℓ / L) atol = 1e-12

    # So the random ring form with that f, read at c̃ → c, IS the ring form.
    S = (c̃ / 3) * log(L * f_conformal(ℓ / L)) + c₁′
    @test check(
        InfiniteRandomnessEntanglementPBC();
        S̄=S,
        c̃=c̃,
        L=L,
        f=f_conformal(ℓ / L),
        c₁′=c₁′,
        atol=1e-12,
    )
    @test check(CFTEntanglementPBC(); S=S, c=c̃, L=L, ℓ=ℓ, c₁=c₁′, atol=1e-12)

    # A second harmonic is what a random chain may carry and a conformal one may
    # not, and it moves the entropy, so the two forms are not interchangeable
    # once it is present.
    A₃ = 0.02
    f_two(v) = (1 - A₃ * 3π) / π * sin(π * v) + A₃ * sin(3π * v)
    @test (1 - A₃ * 3π) / π * π + A₃ * 3π ≈ 1 atol = 1e-14   # normalisation held
    S_two = (c̃ / 3) * log(L * f_two(ℓ / L)) + c₁′
    @test !isapprox(S_two, S; atol=1e-4)
    @test check(
        InfiniteRandomnessEntanglementPBC();
        S̄=S_two,
        c̃=c̃,
        L=L,
        f=f_two(ℓ / L),
        c₁′=c₁′,
        atol=1e-12,
    )
    @test !check(CFTEntanglementPBC(); S=S_two, c=c̃, L=L, ℓ=ℓ, c₁=c₁′, atol=1e-4)
end

@testset "c̃ is not c: the two slopes cannot be read for each other" begin
    rel = InfiniteRandomnessEntanglementSlope()
    @test variable_types(rel) == (EffectiveCentralCharge,)
    @test EffectiveCentralCharge in quantities(rel)

    # The fixed point is not conformal, so an infinite-randomness slope must not
    # produce a CentralCharge row, nor a CFT slope an EffectiveCentralCharge one.
    @test !(CentralCharge in quantities(rel))
    @test !(EffectiveCentralCharge in quantities(CFTEntanglementSlope()))

    # The docstrings promise the entropy arrives through `also_constrains`, so the
    # six new laws must be reachable from the quantity they bound. Without the
    # registration they are silently absent from that lookup and the promise is
    # prose only.
    for r in (
        CFTEntanglementPBC(),
        CFTEntanglementOBC(),
        OffCriticalEntanglementSaturation(),
        HalvedChainEntropyDifference(),
    )
        @test VonNeumannEntropy in AbstractQAtlas.also_constrains(r)
        @test any(x -> x isa typeof(r), relations_constraining(VonNeumannEntropy))
    end

    # The infinite-randomness pair keys on the REDUCTION, following the block in
    # quantity_links.jl: both sources state the disorder average, and a typical
    # sample carries no logarithm, so a bare key would claim the law for the
    # wrong one.
    for r in (InfiniteRandomnessEntanglementSlope(), InfiniteRandomnessEntanglementPBC())
        @test AbstractQAtlas.also_constrains(r) == (DisorderAveraged{VonNeumannEntropy},)
        @test !(VonNeumannEntropy in AbstractQAtlas.also_constrains(r))
    end

    # No dimension slot either: above 1D the entropy is an area law and the
    # fixed point may not even be reached, so there is no family for `d` to
    # index.  A future edit adding one would be claiming a generalisation.
    @test !(SpatialDimension in quantities(rel))

    # The arithmetic is identical, so at the keyword surface nothing separates the
    # two but which keyword the caller types: `name::Type` in `@relation` is
    # bag-key metadata and the generated kernel's kwarg is untyped.
    @test check(rel; dS_dlogℓ=log(2) / 6, c̃=log(2) / 2, ncuts=2, atol=1e-14)
    @test check(
        CFTEntanglementSlope(); dS_dlogℓ=log(2) / 6, c=log(2) / 2, ncuts=2, atol=1e-14
    )

    # The bag is where the separation is a guard rather than a naming convention:
    # the same number keyed as a CentralCharge cannot reach the c̃ slot, so a CFT
    # measurement never picks up the infinite-randomness law, nor the reverse.
    b_cft = bag(CentralCharge => log(2) / 2)
    b_irfp = bag(EffectiveCentralCharge => log(2) / 2)
    slope = (; dS_dlogℓ=log(2) / 6, ncuts=2)
    @test CFTEntanglementSlope() in applicable_relations(b_cft; slope...)
    @test !any(r -> r isa typeof(rel), applicable_relations(b_cft; slope...))
    @test rel in applicable_relations(b_irfp; slope...)
    @test !any(r -> r isa CFTEntanglementSlope, applicable_relations(b_irfp; slope...))
end

function covariance(N)
    A = zeros(2N, 2N)
    for j in 1:N
        A[2j - 1, 2j] = -2.0
    end
    for j in 1:(N - 1)
        A[2j, 2j + 1] = -2.0
    end
    A = A - transpose(A)
    F = svd(A)
    Γ = F.U * F.Vt
    return (Γ .- transpose(Γ)) ./ 2
end
function ed_entropy(Γ, sites)
    idx = vcat(([2j - 1, 2j] for j in sites)...)
    ν = eigvals(Hermitian(im .* Γ[idx, idx]))
    return free_fermion_entanglement_entropy([(1 + real(x)) / 2 for x in ν if real(x) > 0])
end

@testset "the open-chain form against exact diagonalisation and a published constant" begin
    # The closed forms were otherwise checked only by retyping them in the test, so a
    # prefactor transcribed wrongly from the source would be transcribed wrongly here
    # too and nothing would notice. This anchors one of them outside that loop: the
    # entropy comes from an exact free-fermion ground state computed here, and the
    # constant it must reproduce is Iglói & Lin's own measured c₁.
    c = 1 / 2
    c₁_source = log(2) * 0.6904133      # Table 1, converted from their bits to nats

    N = 256
    Γ = covariance(N)
    # Read c₁ back THROUGH the relation, so the production chord is what is exercised.
    recovered(ℓ) =
        solve(CFTEntanglementOBC(), Val(:c₁); S=ed_entropy(Γ, 1:ℓ), c=c, L=N, ℓ=ℓ, ln_g=0.0)

    @test recovered(64) ≈ c₁_source atol = 0.01
    @test recovered(32) ≈ c₁_source atol = 0.01

    # The gap is a finite-size correction, so it has to shrink with ℓ, which a wrong
    # constant would not do.
    @test abs(recovered(64) - c₁_source) < abs(recovered(16) - c₁_source)

    # And the anchor discriminates: dropping the 2 from the open chain's `2L/π` moves
    # the recovered constant by 0.11, twenty times the residual finite-size error, so
    # this test would fail on that transcription rather than absorb it.
    wrong(ℓ) = 2 * (ed_entropy(Γ, 1:ℓ) - (c / 6) * log((N / π) * sin(π * ℓ / N)))
    @test abs(wrong(64) - c₁_source) > 0.1
end

@testset "the chord slope is in domain over the whole chain and ln ℓ is not" begin
    # Exact ED ground state of the critical uniform chain, block at an open end.
    function fitslope(xs, ys)
        x̄, ȳ = sum(xs) / length(xs), sum(ys) / length(ys)
        return sum((xs .- x̄) .* (ys .- ȳ)) / sum((xs .- x̄) .^ 2)
    end
    c, ncuts = 1 / 2, 1

    function slopes(N)
        Γ = covariance(N)
        ls = filter(
            l -> 2 <= l <= N - 2,
            unique(
                round.(
                    Int,
                    N .* [0.03, 0.05, 0.08, 0.12, 0.19, 0.3, 0.45, 0.6, 0.75, 0.88, 0.95],
                ),
            ),
        )
        ys = [ed_entropy(Γ, 1:l) for l in ls]
        chord = [(N / π) * sin(π * l / N) for l in ls]
        idx = findall(l -> l <= N ÷ 4, ls)
        near, yn = ls[idx], ys[idx]
        cn = chord[idx]
        return (
            plain_full=fitslope(log.(ls), ys),
            chord_full=fitslope(log.(chord), ys),
            plain_near=fitslope(log.(near), yn),
            chord_near=fitslope(log.(cn), yn),
        )
    end
    s64, s128, s256 = slopes(64), slopes(128), slopes(256)

    # In the regime the present relation is for, it is right.
    @test isapprox(s128.plain_near, ncuts * c / 6; atol=0.005)

    # Over the whole chain it is not, and by a lot: the block's complement is a few
    # sites, purity caps `S`, and no `ncuts` describes that.
    @test abs(s128.plain_full - ncuts * c / 6) > 0.04

    # The chord form is in domain over the same whole chain.
    @test isapprox(s128.chord_full, ncuts * c / 6; atol=0.01)

    # Two ratios, not one: a 1/L correction halves per doubling and keeps doing it,
    # a wrong law does not move. A threshold cannot separate those.
    e = [abs(s.chord_full - c / 6) for s in (s64, s128, s256)]
    @test e[2] < 0.6 * e[1]
    @test e[3] < 0.6 * e[2]
    @test e[3] < e[2] < e[1]
    q = [abs(s.plain_full - c / 6) for s in (s64, s128, s256)]
    @test q[2] > 0.9 * q[1]
    @test q[3] > 0.9 * q[2]

    # The relation itself: the chord slope reads the central charge back.
    @test isapprox(
        solve(
            CFTEntanglementChordSlope(), Val(:c); dS_dlogchord=s128.chord_full, ncuts=ncuts
        ),
        c;
        atol=0.06,
    )
    # The limit: for ℓ ≪ L the two abscissas coincide, so one dataset gives both
    # relations the same charge.
    @test isapprox(s128.plain_near, s128.chord_near; rtol=0.05)
    @test isapprox(
        solve(CFTEntanglementSlope(), Val(:c); dS_dlogℓ=s128.plain_near, ncuts=ncuts),
        solve(
            CFTEntanglementChordSlope(), Val(:c); dS_dlogchord=s128.chord_near, ncuts=ncuts
        );
        rtol=0.05,
    )

    # With no cuts the residual does not depend on `c`, so both refuse.
    @test_throws "ncuts = 0" residual(
        CFTEntanglementChordSlope(); dS_dlogchord=0.0, c=9.9, ncuts=0
    )
    @test_throws "ncuts = 0" residual(CFTEntanglementSlope(); dS_dlogℓ=0.0, c=9.9, ncuts=0)
end

@testset "every route from an entropy to c returns the same c" begin
    # `related_quantities(CentralCharge)` carries five distinct edges to
    # VonNeumannEntropy, so one measurement can be read for `c` five ways. Each was
    # only ever checked against its own hand-written fixture; here they are made to
    # agree with each other, which no single relation can satisfy alone. A coefficient
    # mistyped in one of them is denied by the other four even where its own test
    # still passes.
    c, c₁, N = 1 / 2, 0.4785, 4096.0
    chord(L, ℓ) = (L / π) * sin(π * ℓ / L)
    ring(ℓ, L=N) = (c / 3) * log(chord(L, ℓ)) + c₁
    open_(ℓ, L=N) = (c / 6) * log(2 * chord(L, ℓ)) + c₁ / 2
    line(ℓ) = (c / 3) * log(ℓ) + c₁

    # A ring: the closed form and the two-length difference.
    c_pbc = solve(CFTEntanglementPBC(), Val(:c); S=ring(1024), L=N, ℓ=1024, c₁=c₁)
    c_halved_ring = solve(
        HalvedChainEntropyDifference(),
        Val(:c);
        ΔS=ring(N / 2) - ring(N / 4, N / 2),
        ncuts=2,
    )
    @test c_pbc ≈ c atol = 1e-12
    @test c_halved_ring ≈ c atol = 1e-12

    # An open chain: same two routes, one cut, and the boundary entropy in the way.
    c_obc = solve(
        CFTEntanglementOBC(), Val(:c); S=open_(1024), L=N, ℓ=1024, c₁=c₁, ln_g=0.0
    )
    c_halved_open = solve(
        HalvedChainEntropyDifference(),
        Val(:c);
        ΔS=open_(N / 2) - open_(N / 4, N / 2),
        ncuts=1,
    )
    @test c_obc ≈ c atol = 1e-12
    @test c_halved_open ≈ c atol = 1e-12

    # An infinite chain: the closed form and the slope through two blocks.
    c_inf = solve(CFTEntanglementInfinite(), Val(:c); S=line(300), ℓ=300, c₁=c₁)
    c_slope = solve(
        CFTEntanglementSlope(),
        Val(:c);
        dS_dlogℓ=(line(300) - line(100)) / (log(300) - log(100)),
        ncuts=2,
    )
    @test c_inf ≈ c atol = 1e-12
    @test c_slope ≈ c atol = 1e-12

    # Off criticality the same c is read from a length that is not the region's.
    ξ = 40.0
    c_sat = solve(
        OffCriticalEntanglementSaturation(), Val(:c); S=2 * (c / 6) * log(ξ), ξ=ξ, ncuts=2
    )
    @test c_sat ≈ c atol = 1e-12

    # All six routes agree to machine precision, which is the claim; the spread a
    # single mistyped coefficient would open is a factor of two, not 1e-12.
    routes = [c_pbc, c_halved_ring, c_obc, c_halved_open, c_inf, c_slope, c_sat]
    @test maximum(routes) - minimum(routes) < 1e-12

    # The spread a single wrong coefficient opens, for scale: the same ring datum
    # read at one cut instead of two returns 2c, so the routes would disagree by c
    # rather than by 1e-12.
    @test solve(
        HalvedChainEntropyDifference(),
        Val(:c);
        ΔS=ring(N / 2) - ring(N / 4, N / 2),
        ncuts=1,
    ) ≈ 2c atol = 1e-12

    # And the ring becomes the infinite chain where its chord does: the two closed
    # forms are one law, not two that happen to agree at one point.
    @test ring(8) ≈ line(8) atol = 2e-5
    @test abs(ring(8) - line(8)) < abs(ring(512) - line(512))
end

@testset "both routes to the effective central charge agree" begin
    # The random side has two edges to EffectiveCentralCharge, and the ring form must
    # reduce to the slope on the same data rather than being a second unrelated law.
    c̃, c₁′, L = log(2) / 2, 0.31, 2048.0
    conf(v) = sin(π * v) / π
    ring̃(ℓ) = (c̃ / 3) * log(L * conf(ℓ / L)) + c₁′
    linẽ(ℓ) = (c̃ / 3) * log(ℓ) + c₁′

    c̃_ring = solve(
        InfiniteRandomnessEntanglementPBC(),
        Val(:c̃);
        S̄=ring̃(512),
        L=L,
        f=conf(512 / L),
        c₁′=c₁′,
    )
    c̃_slope = solve(
        InfiniteRandomnessEntanglementSlope(),
        Val(:c̃);
        dS_dlogℓ=(linẽ(300) - linẽ(100)) / (log(300) - log(100)),
        ncuts=2,
    )
    @test c̃_ring ≈ c̃ atol = 1e-12
    @test c̃_slope ≈ c̃ atol = 1e-12
    @test c̃_ring ≈ c̃_slope atol = 1e-12

    # The two networks must not meet: the same numbers read through the clean route
    # give a different constant, which is the whole reason c̃ is its own quantity.
    @test !isapprox(
        solve(
            CFTEntanglementSlope(),
            Val(:c);
            dS_dlogℓ=(linẽ(300) - linẽ(100)) / (log(300) - log(100)),
            ncuts=2,
        ),
        1 / 2;
        atol=1e-6,
    )
end

@testset "the graph edges are real, and say which traversal can use them" begin
    # Every edge these relations add must name a relation that genuinely links the
    # two quantities, or the graph is advertising a route that is not there.
    edges = related_quantities(CentralCharge)
    named = Set(e.detail for e in edges)
    for r in (
        "CFTEntanglementPBC",
        "CFTEntanglementOBC",
        "CFTEntanglementSlope",
        "HalvedChainEntropyDifference",
        "OffCriticalEntanglementSaturation",
    )
        @test r in named
    end
    @test any(
        e ->
            e.detail == "OffCriticalEntanglementSaturation" &&
            CorrelationLength in (e.from, e.to),
        edges,
    )

    # The edge is a statement about the law, not a promise that `derive` can walk
    # it. The entropy enters as a supplied slot and lives in a bag under a region,
    # so the type-keyed derivation has nothing to match and says so rather than
    # guessing. `CFTEntanglementSlope` behaves the same way, and has since before
    # these relations existed.
    c, c₁, N = 0.5, 0.4785, 4096.0
    b = bag(VonNeumannEntropy => (c / 3) * log((N / π) * sin(π * 1024 / N)) + c₁)
    @test isempty(setdiff(derivable(b; L=N, ℓ=1024.0, c₁=c₁), keys(b)))
    @test_throws "not reachable from the bag" derive(CentralCharge, b; L=N, ℓ=1024.0, c₁=c₁)

    # `finite_size_entropy_report` is the traversal that does reach them, off the
    # region-keyed entry the bag actually holds, and it agrees with the direct
    # solve: the same data passes there and returns c here.
    breg = bag(
        entanglement_entropy(Region(1:1024...)) =>
            (c / 3) * log((N / π) * sin(π * 1024 / N)) + c₁,
        CentralCharge => c,
    )
    rows = finite_size_entropy_report(breg, PBC(Int(N)); c₁=c₁, atol=1e-12)
    @test only(rows).pass
    @test solve(
        CFTEntanglementPBC(),
        Val(:c);
        S=(c / 3) * log((N / π) * sin(π * 1024 / N)) + c₁,
        L=N,
        ℓ=1024,
        c₁=c₁,
    ) ≈ c atol = 1e-12
end

@testset "CFTEntanglementSlope is type-keyed like its cft.jl siblings" begin
    rel = CFTEntanglementSlope()

    # `c` is the typed subject, as in CasimirCentralCharge / CardyDensityOfStates.
    @test variable_types(rel) == (CentralCharge,)
    @test CentralCharge in quantities(rel)
    @test rel in relations_constraining(CentralCharge)

    # Typing `c` does not lose the entropy: it enters through the SUPPLIED
    # derivative, so it has no slot and is declared via `also_constrains`.
    @test also_constrains(rel) == (VonNeumannEntropy,)
    @test VonNeumannEntropy in quantities(rel)
    @test rel in relations_constraining(VonNeumannEntropy)

    # `ncuts` names no quantity, so it stays a supplied coordinate like `L` —
    # it must not have become a bag-visible subject.
    @test :ncuts in variables(rel)
    @test length(variable_types(rel)) == 1
end

@testset "Page average entropy: exact small cases + symmetry" begin
    # two qubits (dA=dB=2): ⟨S⟩ = 1/3 + 1/4 − 1/4 = 1/3, exactly known
    @test page_average_entropy(2, 2) ≈ 1 / 3
    # symmetric in the two dimensions
    @test page_average_entropy(2, 64) == page_average_entropy(64, 2)
    # nearly maximal: ⟨S⟩ → ln(dA) as dB → ∞, with a small positive deficit
    # of order dA/(2 dB) (the harmonic-sum correction is the same order, so
    # this is an asymptotic-form sanity, not a tight match).
    @test page_average_entropy(2, 128) < log(2)                 # below maximal
    @test log(2) - page_average_entropy(2, 128) < 2 / (2 * 128) * 1.5   # deficit ~ dA/(2dB)
    @test page_average_entropy(2, 4096) ≈ log(2) rtol = 1e-3     # deficit vanishes as dB→∞
    @test_throws ErrorException page_average_entropy(0, 4)
end

@testset "Page formula vs a Haar-random-state average" begin
    # dA=2, dB=3: sample random pure states, average the reduced-ρ_A entropy
    dA, dB = 2, 3
    rng = MersenneTwister(2024)
    Ssum = 0.0
    M = 4000
    for _ in 1:M
        ψ = randn(rng, ComplexF64, dA * dB)
        ψ ./= norm(ψ)
        Ψ = reshape(ψ, dA, dB)              # A ⊗ B
        ρA = Ψ * Ψ'                         # reduced density matrix on A
        ev = filter(>(1e-14), real(eigvals(Hermitian(ρA))))
        Ssum += -sum(ev .* log.(ev))
    end
    @test isapprox(Ssum / M, page_average_entropy(dA, dB); rtol=0.03)
end

@testset "entropy zoo: Rényi/Tsallis moments, mutual & conditional, Klein" begin
    using AbstractQAtlas: check, solve, slack, residual, AbstractInequality

    # Schmidt spectrum {p, 1−p}: moments Tr ρ^α = p^α + (1−p)^α
    p = 0.3
    mom(α) = p^α + (1 - p)^α
    SvN = -p * log(p) - (1 - p) * log(1 - p)

    # Rényi from the moment reduces to −ln(purity) at α=2 (= RenyiTwoPurity) …
    S2 = solve(RenyiEntropyMoment(), Val(:Sα); moment=mom(2), α=2)
    @test S2 ≈ -log(mom(2)) atol = 1e-12
    @test check(RenyiTwoPurity(); S2=S2, purity=mom(2), atol=1e-12)      # cross-relation
    # … and → S_vN as α → 1 (l'Hôpital limit, numerically)
    @test solve(RenyiEntropyMoment(), Val(:Sα); moment=mom(1.0001), α=1.0001) ≈ SvN atol =
        1e-3

    # Tsallis from the moment; q → 1 limit is S_vN
    @test check(
        TsallisEntropyMoment(); Sq=(1 - mom(2)) / (2 - 1), moment=mom(2), q=2, atol=1e-12
    )
    @test solve(TsallisEntropyMoment(), Val(:Sq); moment=mom(1.0001), q=1.0001) ≈ SvN atol =
        1e-3

    # mutual information I = S_A+S_B−S_AB ≥ 0, and on a pure state (S_AB=0)
    # with S_A=S_B=H(p): I = 2 H(p); consistency with Subadditivity slack
    S_A = S_B = SvN
    S_AB = 0.0
    I = solve(MutualInformationDefinition(), Val(:I); S_A=S_A, S_B=S_B, S_AB=S_AB)
    @test I ≈ 2 * SvN atol = 1e-12
    @test slack(Subadditivity(); S_A=S_A, S_B=S_B, S_AB=S_AB) ≈ I atol = 1e-12   # I = subadditivity slack

    # conditional entropy of a pure entangled state is NEGATIVE: S(A|B)=S_AB−S_B=−S_B
    S_cond = solve(ConditionalEntropyDefinition(), Val(:S_cond); S_AB=S_AB, S_B=S_B)
    @test S_cond ≈ -S_B atol = 1e-12
    @test S_cond < 0                                                     # entanglement witness

    # Klein's inequality S(ρ‖σ) ≥ 0 (inequality kind), zero iff ρ=σ
    @test RelativeEntropyNonNegativity() isa AbstractInequality
    @test check(RelativeEntropyNonNegativity(); S_rel=0.4)
    @test slack(RelativeEntropyNonNegativity(); S_rel=0.0) == 0.0        # ρ = σ saturates
    @test !check(RelativeEntropyNonNegativity(); S_rel=-1e-3, atol=1e-9)
end

@testset "measurement + Markov entropies on concrete states" begin
    using AbstractQAtlas: check, solve, slack, AbstractInequality
    using LinearAlgebra: eigvals, diagm, Hermitian

    # a 2×2 density matrix with coherences; measure (dephase) in the z basis
    a, c = 0.7, 0.3               # populations
    off = 0.35                    # coherence (|off|² ≤ a·c for positivity: 0.1225 ≤ 0.21 ✓)
    ρ = [a off; off c]
    Δρ = [a 0.0; 0.0 c]           # dephased (diagonal part)
    ent(M) = (λ=filter(>(1e-15), real(eigvals(Hermitian(M)))); -sum(x -> x * log(x), λ))
    S = ent(ρ)
    S_meas = ent(Δρ)
    # relative entropy S(ρ‖Δρ) = Tr ρ(ln ρ − ln Δρ)
    using LinearAlgebra: tr
    lnρ = log(Hermitian(ρ))
    lnΔρ = diagm([log(a), log(c)])
    S_rel = real(tr(ρ * (lnρ - lnΔρ)))

    # measurement does not decrease entropy: S(Δρ) ≥ S(ρ)
    @test check(MeasurementEntropyIncrease(); S_meas=S_meas, S=S)
    @test slack(MeasurementEntropyIncrease(); S_meas=S_meas, S=S) > 0     # strict, coherences present
    # …and the gain equals the relative entropy to the dephased state (exact identity)
    @test check(MeasurementEntropyRelative(); S_meas=S_meas, S=S, S_rel=S_rel, atol=1e-10)
    @test (S_meas - S) ≈ S_rel atol = 1e-10
    # a state already diagonal saturates the inequality (Δρ = ρ)
    @test slack(MeasurementEntropyIncrease(); S_meas=S, S=S) == 0.0

    # MarkovEntropy = conditional mutual information = strong-subadditivity slack
    S_A, S_B, S_C = 0.4, 0.7, 0.5
    S_AB, S_BC, S_ABC = S_A + S_B, S_B + S_C, S_A + S_B + S_C          # product (Markov) state
    I_cmi = solve(
        MarkovEntropyDefinition(), Val(:I_cmi); S_AB=S_AB, S_BC=S_BC, S_ABC=S_ABC, S_B=S_B
    )
    @test I_cmi ≈ slack(StrongSubadditivity(); S_AB=S_AB, S_BC=S_BC, S_ABC=S_ABC, S_B=S_B) atol =
        1e-12
    @test I_cmi ≈ 0 atol = 1e-12                                        # a Markov chain: I(A:C|B)=0

    @test MeasurementEntropyIncrease() isa AbstractInequality
    @test tensor_rank(MeasurementEntropy()) == 0
    @test tensor_rank(MarkovEntropy()) == 0
end

@testset "free-fermion entanglement from the correlation matrix (Peschel)" begin
    using AbstractQAtlas:
        check, solve, free_fermion_entanglement_entropy, free_fermion_renyi_entropy
    using LinearAlgebra: eigvals, Hermitian

    # GENUINE Gaussian case: the bonding orbital (c†₁+c†₂)/√2|0⟩ of a 2-site
    # hopping model.  Correlation matrix C_ij = ⟨c†_i c_j⟩ = 1/2 (all entries);
    # trace out site 2 ⇒ C_A = [1/2], eigenvalue ζ = 1/2 ⇒ S_A = ln 2 (one Bell pair).
    C = [0.5 0.5; 0.5 0.5]
    C_A = C[1:1, 1:1]                              # region A = site 1
    ζ = real(eigvals(Hermitian(C_A)))
    @test ζ ≈ [0.5] atol = 1e-12
    @test free_fermion_entanglement_entropy(ζ) ≈ log(2) atol = 1e-12

    # a filled/empty mode contributes nothing; a maximal mode ln 2
    @test free_fermion_entanglement_entropy([0.0, 1.0, 1.0]) == 0.0
    @test free_fermion_entanglement_entropy([0.5, 0.5, 0.5]) ≈ 3 * log(2) atol = 1e-12

    # Peschel single-particle spectrum ε = ln((1−ζ)/ζ): ζ=1/2 ⇒ ε=0 (max entangled)
    @test check(EntanglementSpectrumCorrelation(); ε=0.0, ζ=0.5, atol=1e-12)
    @test solve(EntanglementSpectrumCorrelation(), Val(:ε); ζ=0.3) ≈ log(0.7 / 0.3)
    # inverting the spectrum recovers the Fermi-Dirac occupation ζ = 1/(e^ε+1)
    ε = 1.1
    ζinv = 1 / (exp(ε) + 1)
    @test check(EntanglementSpectrumCorrelation(); ε=ε, ζ=ζinv, atol=1e-12)

    # Rényi → von Neumann as n → 1 (INDEPENDENT: two different formulas agree)
    ζset = [0.15, 0.5, 0.82, 0.97]
    @test free_fermion_renyi_entropy(ζset, 1.0001) ≈ free_fermion_entanglement_entropy(ζset) atol =
        1e-3
    # Rényi-2 = −Σ ln(ζ²+(1−ζ)²) matches the moment definition per mode
    @test free_fermion_renyi_entropy([0.3], 2) ≈ -log(0.3^2 + 0.7^2) atol = 1e-12
    @test_throws ErrorException free_fermion_renyi_entropy(ζset, 1)   # n=1 is the vN limit
end

@testset "multipartite: monogamy / three-tangle on GHZ and W states" begin
    using AbstractQAtlas: check, solve, slack, AbstractInequality

    # tangle = concurrence²; a Bell pair has C=1 ⇒ τ=1
    @test check(ConcurrenceTangle(); τ=1.0, C=1.0, atol=1e-12)
    @test solve(ConcurrenceTangle(), Val(:τ); C=2 / 3) ≈ 4 / 9   # W-state pair concurrence

    # GHZ = (|000⟩+|111⟩)/√2: pairwise reduced states are separable ⇒ τ_AB=τ_AC=0,
    # while A is maximally entangled with BC ⇒ τ(A:BC)=1, so the three-tangle = 1.
    τ_ABC_ghz, τ_AB_ghz, τ_AC_ghz = 1.0, 0.0, 0.0
    @test check(Monogamy(); τ_ABC=τ_ABC_ghz, τ_AB=τ_AB_ghz, τ_AC=τ_AC_ghz)
    τ3_ghz = solve(
        ThreeTangleDefinition(), Val(:τ3); τ_ABC=τ_ABC_ghz, τ_AB=τ_AB_ghz, τ_AC=τ_AC_ghz
    )
    @test τ3_ghz ≈ 1.0 atol = 1e-12                              # GHZ: genuine tripartite entanglement
    @test τ3_ghz ≈ slack(Monogamy(); τ_ABC=τ_ABC_ghz, τ_AB=τ_AB_ghz, τ_AC=τ_AC_ghz) atol =
        1e-12

    # W = (|001⟩+|010⟩+|100⟩)/√3: τ(A:B)=τ(A:C)=4/9, τ(A:BC)=8/9 ⇒ three-tangle = 0
    # (W saturates monogamy — no residual tripartite tangle).
    τ_ABC_w, τ_AB_w, τ_AC_w = 8 / 9, 4 / 9, 4 / 9
    @test check(Monogamy(); τ_ABC=τ_ABC_w, τ_AB=τ_AB_w, τ_AC=τ_AC_w, atol=1e-12)
    @test slack(Monogamy(); τ_ABC=τ_ABC_w, τ_AB=τ_AB_w, τ_AC=τ_AC_w) ≈ 0 atol = 1e-12
    @test solve(
        ThreeTangleDefinition(), Val(:τ3); τ_ABC=τ_ABC_w, τ_AB=τ_AB_w, τ_AC=τ_AC_w
    ) ≈ 0 atol = 1e-12

    # monogamy is a genuine bound: over-sharing (τ_AB+τ_AC > τ_ABC) is forbidden
    @test Monogamy() isa AbstractInequality
    @test !check(Monogamy(); τ_ABC=0.5, τ_AB=0.4, τ_AC=0.4, atol=1e-9)
end

@testset "tripartite information + Kitaev–Preskill TEE" begin
    using AbstractQAtlas: check, solve
    # I₃ = I(A:B) + I(A:C) − I(A:BC)
    @test check(
        TripartiteInformationDefinition();
        I3=0.3 + 0.5 - 1.1,
        I_AB=0.3,
        I_AC=0.5,
        I_ABC=1.1,
        atol=1e-12,
    )
    @test solve(
        TripartiteInformationDefinition(), Val(:I3); I_AB=0.3, I_AC=0.5, I_ABC=1.1
    ) ≈ -0.3

    # Kitaev–Preskill: the alternating tripartite sum isolates −γ. Build a
    # toric-code-like assignment (all singles a, pairs b, triple c):
    # 3a − 3b + c = −γ ⇒ pick a,b,c so γ = ln 2 (toric code total dimension D=2).
    a, b = 1.7, 2.4
    γ = log(2)
    c = -γ - (3a - 3b)                          # so 3a − 3b + c = −γ
    @test check(
        KitaevPreskillTEE();
        γ=γ,
        S_A=a,
        S_B=a,
        S_C=a,
        S_AB=b,
        S_BC=b,
        S_CA=b,
        S_ABC=c,
        atol=1e-12,
    )
    @test solve(
        KitaevPreskillTEE(), Val(:γ); S_A=a, S_B=a, S_C=a, S_AB=b, S_BC=b, S_CA=b, S_ABC=c
    ) ≈ log(2) atol = 1e-12                      # extracts γ = ln 2

    @test tensor_rank(Concurrence()) == 0
    @test tensor_rank(TopologicalEntanglementEntropy()) == 0
end
