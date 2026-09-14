# test/relations/test_region_entropy.jl — §5 Region set layer + region-keyed
# entanglement-entropy auto-discovery (design doc §5/§8b, Phase-2 P1).

using AbstractQAtlas
using Test

struct _NoFormBC <: AbstractQAtlas.BoundaryCondition end
const AQ = AbstractQAtlas

@testset "Region set algebra (dimension-agnostic)" begin
    A, B = Region(1, 2), Region(3, 4)
    @test disjoint(A, B)
    @test !disjoint(A, Region(2, 3))
    @test A ∪ B == Region(1, 2, 3, 4)
    @test A ∩ Region(2, 3) == Region(2)
    @test Region(1) ⊆ A && !(Region(1, 5) ⊆ A)
    @test length(A) == 2 && !isempty(A) && isempty(Region{Int}(Set{Int}()))
    # value-based equality/hash — order-insensitive
    @test Region(2, 1) == Region(1, 2)
    @test hash(Region(2, 1)) == hash(Region(1, 2))
    # ND sites: 2D tuples work identically (the set layer is ND from day one)
    a2, b2 = Region((1, 1), (1, 2)), Region((2, 1))
    @test disjoint(a2, b2)
    @test a2 ∪ b2 == Region((1, 1), (1, 2), (2, 1))
end

@testset "entanglement_entropy() region-keyed bag" begin
    @test entanglement_entropy(1, 2) == entanglement_entropy(Region(1, 2))                 # value-based key
    @test entanglement_entropy(1, 2) isa VariableKey
    @test entanglement_entropy(1, 2).type === VonNeumannEntropy
    @test entanglement_entropy(1, 2).support isa RegionSupport
    @test entanglement_entropy(2, 1) == entanglement_entropy(1, 2)                          # region order-insensitive
    @test bag(entanglement_entropy(1) => 0.5, entanglement_entropy(1, 2) => 1.0)[entanglement_entropy(
        1
    )] == 0.5
end

@testset "region_report: auto-discovery of subadditivity + Araki–Lieb" begin
    # subadditive bag: S(A)+S(B) ≥ S(A∪B) and S(A∪B) ≥ |S(A)−S(B)| — all hold
    good = bag(
        entanglement_entropy(1) => 0.7,
        entanglement_entropy(2) => 0.7,
        entanglement_entropy(1, 2) => 1.0,
    )
    rep = region_report(good)
    @test length(rep) == 2                                       # Subadditivity + ArakiLieb, one region pair
    @test all(r -> r.pass, rep)
    @test region_check_all(good)
    @test Set(nameof(typeof(r.relation)) for r in rep) == Set((:Subadditivity, :ArakiLieb))

    # a NEGATIVE mutual information (unphysical — a broken calc) is caught, on the
    # right region pair, with no A/B/AB hand-labeling
    bad = bag(
        entanglement_entropy(1) => 0.5,
        entanglement_entropy(2) => 0.5,
        entanglement_entropy(1, 2) => 1.5,
    )   # I(A:B) = −0.5
    @test !region_check_all(bad)
    viol = [r for r in region_report(bad) if !r.pass]
    @test length(viol) == 1
    @test viol[1].relation isa Subadditivity
    @test Set(viol[1].regions) == Set((Region(1), Region(2)))
    @test viol[1].slack ≈ -0.5

    # 3 regions ⇒ every disjoint pair whose union is present is auto-discovered
    b3 = bag(
        entanglement_entropy(1) => 0.6,
        entanglement_entropy(2) => 0.6,
        entanglement_entropy(3) => 0.6,
        entanglement_entropy(1, 2) => 1.0,
        entanglement_entropy(1, 3) => 1.0,
        entanglement_entropy(2, 3) => 1.0,
    )
    # 3 disjoint pairs × 2 pair-relations + 3 weak-monotonicity triples (one per middle;
    # no full-system entropy so no SSA)
    @test length(region_report(b3)) == 9
    @test region_check_all(b3)

    # empty match ⇒ false, never a silent green
    @test !region_check_all(bag(entanglement_entropy(1) => 0.5))
    @test isempty(region_report(bag(entanglement_entropy(1) => 0.5)))
    # a non-disjoint pair (or a missing union / missing S(B)) yields no instance
    @test isempty(
        region_report(
            bag(entanglement_entropy(1) => 0.5, entanglement_entropy(1, 2) => 1.0)
        ),
    )
end

@testset "region_report: strong subadditivity (triples) + mutual_information" begin
    ee = entanglement_entropy
    # SSA-valid bag (atoms A={1}, B={2}, C={3}): S(A∪B)+S(B∪C) ≥ S(A∪B∪C)+S(B)
    b = bag(
        ee(1) => 0.5,
        ee(2) => 0.5,
        ee(3) => 0.5,
        ee(1, 2) => 1.0,
        ee(2, 3) => 1.0,
        ee(1, 2, 3) => 1.2,
    )
    ssa = [r for r in region_report(b) if r.relation isa StrongSubadditivity]
    @test length(ssa) == 1                                # one triple, B={2} the shared middle
    @test ssa[1].pass && ssa[1].slack ≈ 0.3
    @test Set(ssa[1].regions) == Set((Region(1), Region(2), Region(3)))
    @test region_check_all(b)                             # SSA + subadditivity + Araki–Lieb all hold

    # a FULLY-connected triple: all 3 choices of the shared middle B are valid, so the
    # enumeration must emit exactly 3 distinct SSA instances (no dup, no omission)
    full = bag(
        ee(1) => 0.5,
        ee(2) => 0.5,
        ee(3) => 0.5,
        ee(1, 2) => 1.0,
        ee(1, 3) => 1.0,
        ee(2, 3) => 1.0,
        ee(1, 2, 3) => 1.2,
    )
    ssa_full = [r for r in region_report(full) if r.relation isa StrongSubadditivity]
    @test length(ssa_full) == 3                           # one per choice of the middle B
    @test length(unique(r.regions[2] for r in ssa_full)) == 3      # 3 distinct middles
    @test all(r -> r.pass, ssa_full)

    # a broken calc — negative conditional mutual information I(A:C|B) < 0 — is caught
    bad = bag(
        ee(1) => 0.5,
        ee(2) => 0.5,
        ee(3) => 0.5,
        ee(1, 2) => 1.0,
        ee(2, 3) => 1.0,
        ee(1, 2, 3) => 1.8,
    )
    @test !region_check_all(bad)
    viol = [r for r in region_report(bad) if r.relation isa StrongSubadditivity && !r.pass]
    @test length(viol) == 1 && viol[1].slack ≈ -0.3

    # mutual_information helper: I(A:B) = S(A) + S(B) − S(A∪B)
    @test mutual_information(b, Region(1), Region(2)) ≈ 0.0                 # 0.5 + 0.5 − 1.0
    @test mutual_information(
        bag(ee(1) => 0.7, ee(2) => 0.7, ee(1, 2) => 1.0), Region(1), Region(2)
    ) ≈ 0.4
    @test_throws ErrorException mutual_information(bag(ee(1) => 0.5), Region(1), Region(2))
end

@testset "region_report: weak monotonicity (partial-data triples)" begin
    ee = entanglement_entropy
    # weak monotonicity S(A∪B)+S(B∪C) ≥ S(A)+S(C) needs NO full-system S(A∪B∪C), so it is
    # auto-discovered on a triple where SSA cannot even be posed — and catches a broken
    # calc there (a negative WM slack is unphysical, like a negative mutual information)
    wmbad = bag(
        ee(1) => 1.0,
        ee(2) => 0.5,
        ee(3) => 1.0,
        ee(1, 2) => 0.8,
        ee(2, 3) => 0.8,       # middle B={2}: 0.8+0.8 − S(1) − S(3) = −0.4 < 0 (unphysical)
    )
    rep = region_report(wmbad)
    # NO strong-subadditivity instance exists here (S(A∪B∪C) absent) — this is exactly the
    # extra coverage weak monotonicity provides
    @test isempty([r for r in rep if r.relation isa StrongSubadditivity])
    viol = [r for r in rep if !r.pass]
    @test length(viol) == 1
    @test viol[1].relation isa WeakMonotonicity
    @test Set(viol[1].regions) == Set((Region(1), Region(2), Region(3)))
    @test viol[1].slack ≈ -0.4
    @test !region_check_all(wmbad)

    # on a valid full triple, WM and SSA are BOTH discovered and BOTH pass
    good = bag(
        ee(1) => 0.5,
        ee(2) => 0.5,
        ee(3) => 0.5,
        ee(1, 2) => 1.0,
        ee(2, 3) => 1.0,
        ee(1, 2, 3) => 1.2,
    )
    kinds = Set(nameof(typeof(r.relation)) for r in region_report(good))
    @test :WeakMonotonicity in kinds && :StrongSubadditivity in kinds
    @test region_check_all(good)
    wm = only(r for r in region_report(good) if r.relation isa WeakMonotonicity)
    @test wm.slack ≈ 1.0                      # 1.0 + 1.0 − 0.5 − 0.5, middle B={2}
end

@testset "multipartite region helpers: CMI, tripartite info, KP topological EE" begin
    ee = entanglement_entropy
    A, B, C = Region(1), Region(2), Region(3)
    # pairwise-DISTINCT singles and pairs so a positional swap among {S_A,S_B,S_C} or
    # {S_AB,S_AC,S_BC} — the failure mode of the generator destructuring — cannot hide
    b = bag(
        ee(1) => 0.5,
        ee(2) => 0.6,
        ee(3) => 0.7,
        ee(1, 2) => 1.0,
        ee(1, 3) => 1.1,
        ee(2, 3) => 1.2,
        ee(1, 2, 3) => 1.2,
    )
    # I(A:C|B) = S(A∪B)+S(B∪C)−S(A∪B∪C)−S(B) = 1.0+1.2−1.2−0.6   (the SSA slack; S(B)=S₂
    # is distinct from S₁/S₃, so picking the wrong middle region would change the answer)
    cmi = conditional_mutual_information(b, A, B, C)
    @test cmi ≈ 0.4
    # I₃ = S_A+S_B+S_C−S_AB−S_AC−S_BC+S_ABC = 1.8−3.3+1.2 ; Kitaev–Preskill γ = −I₃
    @test tripartite_information(b, A, B, C) ≈ -0.3
    @test topological_entanglement_entropy(b, A, B, C) ≈ 0.3
    @test topological_entanglement_entropy(b, A, B, C) ≈ -tripartite_information(b, A, B, C)
    # CMI agrees with the auto-discovered StrongSubadditivity slack for the same triple
    @test any(r -> r.relation isa StrongSubadditivity && r.slack ≈ cmi, region_report(b))
    # a missing entropy is a loud error, never silently wrong
    @test_throws ErrorException tripartite_information(bag(ee(1) => 0.5), A, B, C)
    @test_throws ErrorException conditional_mutual_information(
        bag(ee(1, 2) => 1.0), A, B, C
    )
    @test_throws ErrorException topological_entanglement_entropy(bag(ee(1) => 0.5), A, B, C)
    # overlapping regions are NOT a tripartition (unions collapse ⇒ a physical-looking
    # but meaningless number) — caught by the disjointness guard, on a bag where every
    # needed entropy IS present so only the guard can be firing
    @test_throws ErrorException conditional_mutual_information(b, A, B, A)   # C == A
    @test_throws ErrorException tripartite_information(b, A, B, A)
    @test_throws ErrorException topological_entanglement_entropy(b, A, B, A)
end

@testset "region_tee_report: auto-discovered tripartite info + KP topological EE" begin
    ee = entanglement_entropy
    γ = log(2)
    # an area-law-canceling tripartition (as in a Kitaev–Preskill disk split into three
    # sectors): the six area terms cancel in the alternating sum and S(A∪B∪C) is set so
    # ΣₐₗₜS = −γ, so the report isolates γ = ln2 (toric code)
    b = bag(
        ee(1) => 1.0,
        ee(2) => 1.0,
        ee(3) => 1.0,
        ee(1, 2) => 1.5,
        ee(1, 3) => 1.5,
        ee(2, 3) => 1.5,
        ee(1, 2, 3) => 1.5 - γ,     # = 3·1.5 − 3·1.0 − γ ⇒ alternating sum = −γ
    )
    rep = region_tee_report(b)
    @test length(rep) == 1                    # one unordered triple (I₃ symmetric ⇒ no dup)
    @test rep[1].topological_entanglement_entropy ≈ γ          # ln 2
    @test rep[1].tripartite_information ≈ -γ
    @test Set(rep[1].regions) == Set((Region(1), Region(2), Region(3)))
    # the auto-discovered value equals the manual helpers on the same triple (no drift)
    @test rep[1].tripartite_information ≈
        tripartite_information(b, Region(1), Region(2), Region(3))
    @test rep[1].topological_entanglement_entropy ≈
        topological_entanglement_entropy(b, Region(1), Region(2), Region(3))

    # a trivial (product) state has γ = 0: S is strictly additive over the tripartition
    triv = bag(
        ee(1) => 0.4,
        ee(2) => 0.5,
        ee(3) => 0.6,
        ee(1, 2) => 0.9,       # S additive: 0.4+0.5
        ee(1, 3) => 1.0,       # 0.4+0.6
        ee(2, 3) => 1.1,       # 0.5+0.6
        ee(1, 2, 3) => 1.5,    # 0.4+0.5+0.6
    )
    @test only(region_tee_report(triv)).topological_entanglement_entropy ≈ 0.0 atol = 1e-12

    # MULTIPLE simultaneous tripartitions: 4 disjoint atoms with every pair+triple union
    # present ⇒ C(4,3)=4 unordered triples, each discovered exactly once — the multipartite
    # analog of the SSA multiplicity test above (guards the i<j<k enumeration against
    # dup/omission, which the single-triple cases cannot)
    s = Dict(1 => 0.4, 2 => 0.5, 3 => 0.6, 4 => 0.7)
    prod4 = bag(
        (ee(i) => s[i] for i in 1:4)...,
        (
            ee(p...) => sum(s[i] for i in p) for
            p in ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
        )...,
        (
            ee(t...) => sum(s[i] for i in t) for
            t in ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4))
        )...,
    )
    rep4 = region_tee_report(prod4)
    @test length(rep4) == 4                                  # C(4,3): one row per unordered triple
    @test Set(Set(r.regions) for r in rep4) ==
        Set(Set(Region.(t)) for t in ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)))
    @test all(r -> isapprox(r.topological_entanglement_entropy, 0.0; atol=1e-12), rep4)

    # no valid tripartition ⇒ empty, never a silent zero: fewer than 3 disjoint regions…
    @test isempty(region_tee_report(bag(ee(1) => 0.5, ee(2) => 0.5, ee(1, 2) => 1.0)))
    # …or the full-triple entropy S(A∪B∪C) is absent
    nofull = bag(
        ee(1) => 1.0,
        ee(2) => 1.0,
        ee(3) => 1.0,
        ee(1, 2) => 1.5,
        ee(1, 3) => 1.5,
        ee(2, 3) => 1.5,
    )
    @test isempty(region_tee_report(nofull))
end

@testset "Fermionic entanglement entropy is a SEPARATE family" begin
    ee = entanglement_entropy
    fe = fermionic_entanglement_entropy
    A, B = Region(1, 2), Region(5, 6)

    # different key on the same region — the whole point of the split
    @test ee(A) != fe(A)
    @test fe(A).type === FermionicEntanglementEntropy
    @test fe(A).support == RegionSupport(A)
    @test fe(1, 2) == fe(Region(1, 2))

    # opt-in trait: entanglement measures that are NOT strongly subadditive
    # away from the von Neumann limit stay out
    @test AQ.obeys_entropy_inequalities(VonNeumannEntropy)
    @test AQ.obeys_entropy_inequalities(FermionicEntanglementEntropy)
    @test !AQ.obeys_entropy_inequalities(RenyiEntropy)
    @test !AQ.obeys_entropy_inequalities(TsallisEntropy)
    @test !AQ.obeys_entropy_inequalities(MutualInformation)

    # Measured on the open XXZ(Δ = 0) chain at N = 12 (QAtlas, dense ED for the
    # spin entropy and the free-fermion covariance for the fermionic one).  The
    # two agree exactly on the contiguous blocks and differ only on the
    # disconnected union — which is the entire content of the Jordan–Wigner
    # caveat, and the reason the union entropy alone cannot be shared.
    S_spin = Dict(A => 0.519867, B => 0.814639, A ∪ B => 1.112324)
    S_ferm = Dict(A => 0.519867, B => 0.814639, A ∪ B => 1.224109)
    @test S_spin[A] == S_ferm[A] && S_spin[B] == S_ferm[B]
    @test S_spin[A ∪ B] != S_ferm[A ∪ B]

    both = bag((ee(R) => v for (R, v) in S_spin)..., (fe(R) => v for (R, v) in S_ferm)...)
    rep = region_report(both)
    @test !isempty(rep)
    @test all(r -> r.pass, rep)
    @test region_check_all(both)

    # The detector for a merged sweep.  `_region_entropies` returns Region =>
    # value, so if the two families shared one Dict their equal region keys would
    # COLLIDE and one union entropy would silently overwrite the other: the report
    # would come back at single-family length with one family's numbers gone.
    single = region_report(bag((ee(R) => v for (R, v) in S_spin)...))
    @test length(rep) == 2 * length(single)

    # …and both families' own numbers survive, with DIFFERENT slacks on the same
    # region pair — the sharpest statement that nothing was overwritten
    subs = [r.slack for r in rep if r.relation isa Subadditivity && r.regions == (A, B)]
    @test length(subs) == 2
    @test !isapprox(subs[1], subs[2]; atol=1e-6)
    @test Set(round.(subs; digits=6)) == Set(
        round.([0.519867 + 0.814639 - 1.112324, 0.519867 + 0.814639 - 1.224109]; digits=6)
    )

    # `quantity` selects the family for the mutual information
    I_spin = mutual_information(both, A, B)
    I_ferm = mutual_information(both, A, B; quantity=FermionicEntanglementEntropy)
    @test I_spin ≈ 0.222182 atol = 1e-6
    @test I_ferm ≈ 0.110397 atol = 1e-6
    @test I_spin > I_ferm            # the JW string ADDS spin correlation across the gap

    # a fermionic-only bag is swept on its own
    fonly = bag((fe(R) => v for (R, v) in S_ferm)...)
    @test length(region_report(fonly)) == length(single)
    @test region_check_all(fonly)

    # a family that has not opted in is invisible rather than swept — even when
    # its values would violate subadditivity outright
    ronly = bag(
        VariableKey(RenyiEntropy, RegionSupport(A)) => 0.5,
        VariableKey(RenyiEntropy, RegionSupport(B)) => 0.5,
        VariableKey(RenyiEntropy, RegionSupport(A ∪ B)) => 3.0,
    )
    @test isempty(region_report(ronly))
    @test !region_check_all(ronly)   # an empty match is never a silent green

    # row order is deterministic (regions come out of a Dict)
    @test [(typeof(r.relation), r.regions) for r in region_report(both)] == [(typeof(r.relation), r.regions) for r in region_report(both)]
end

@testset "MaxEntropyBound is discovered once the local dimension is supplied" begin
    # Without `local_dim` nothing knows the Hilbert space, so a flatly impossible
    # entropy produces NO row — the gap this opt-in closes.
    impossible = bag(entanglement_entropy(1) => 5.0)     # 5 nats on one qubit
    @test isempty(region_report(impossible))
    @test !region_check_all(impossible; local_dim=2)

    rows = region_report(impossible; local_dim=2)
    @test length(rows) == 1
    @test only(rows).relation isa MaxEntropyBound
    @test only(rows).regions == (Region(1),)
    @test !only(rows).pass
    @test only(rows).slack ≈ log(2) - 5.0                # ln d − S, the deficit

    # A legal bag passes, and the bound is per-region: two sites allow 2 ln 2.
    ok = bag(entanglement_entropy(1) => 0.4, entanglement_entropy(1, 2) => 1.2)
    @test all(r -> r.pass, region_report(ok; local_dim=2))
    @test region_check_all(ok; local_dim=2)
    # ...and it CAN disagree at that size: 1.2 < 2ln2 = 1.386 passes, 1.5 does not.
    @test !region_check_all(
        bag(entanglement_entropy(1) => 0.4, entanglement_entropy(1, 2) => 1.5); local_dim=2
    )

    # Qutrits give a strictly looser bound on the same data, so the dimension is
    # doing work rather than being decoration.
    tight = bag(entanglement_entropy(1) => 1.0)          # ln2 = 0.693 < 1 < ln3 = 1.099
    @test !region_check_all(tight; local_dim=2)
    @test region_check_all(tight; local_dim=3)

    @test_throws ArgumentError region_report(impossible; local_dim=1)

    # Opt-in only: the existing inequality rows are unchanged when it is omitted.
    b = bag(
        entanglement_entropy(1) => 0.7,
        entanglement_entropy(2) => 0.7,
        entanglement_entropy(1, 2) => 1.0,
    )
    @test length(region_report(b; local_dim=2)) == length(region_report(b)) + 3
end

@testset "a pure global state needs no complementarity relation" begin
    # S(A∪B) = 0 makes A∪B pure, so S(A) must equal S(B).  Araki–Lieb already
    # says so — `S_AB ≥ |S_A − S_B|` collapses to equality at S_AB = 0 — which is
    # why there is no separate relation for it.
    violating = bag(
        entanglement_entropy(1) => 0.7,
        entanglement_entropy(2) => 0.3,
        entanglement_entropy(1, 2) => 0.0,
    )
    rows = region_report(violating)
    al = only(filter(r -> r.relation isa ArakiLieb, rows))
    @test !al.pass
    @test al.slack ≈ -0.4                                 # 0 − |0.7 − 0.3|
    @test !region_check_all(violating)

    # The same bag with complementarity RESPECTED passes every row.
    respecting = bag(
        entanglement_entropy(1) => 0.5,
        entanglement_entropy(2) => 0.5,
        entanglement_entropy(1, 2) => 0.0,
    )
    @test region_check_all(respecting)
end

@testset "local_dim: validated at the boundary, and what a pass without it means" begin
    # Validated in `region_report` itself, so an invalid value is rejected even
    # when the bag holds no entropy family to sweep — otherwise a typo'd
    # `local_dim` would look exactly like "no data yet".
    @test_throws ArgumentError region_report(bag(); local_dim=1)
    @test_throws ArgumentError region_report(
        bag(entanglement_entropy(1) => 0.5); local_dim=1
    )
    # Non-integer and non-finite dimensions cannot reach the bound at all: `2.5`
    # would loosen it silently and `Inf` would make `ln d` infinite, passing ANY
    # entropy.  The `Int` annotation refuses both at the call boundary.
    @test_throws TypeError region_report(bag(); local_dim=2.5)
    @test_throws TypeError region_report(bag(); local_dim=Inf)

    # A pass WITHOUT `local_dim` does not include the maximum-entropy bound: the
    # inequalities that were checked can all hold on data the omitted one refuses.
    impossible = bag(
        entanglement_entropy(1) => 5.0,       # 5 nats on one qubit
        entanglement_entropy(2) => 5.0,
        entanglement_entropy(1, 2) => 8.0,
    )
    @test region_check_all(impossible)                     # subadditive, Araki–Lieb: true
    @test !region_check_all(impossible; local_dim=2)       # ...and impossible
end

@testset "finite-size forms are reached from a region-keyed bag" begin
    c, c₁, N = 0.5, 0.4785, 64
    ring(ℓ) = (c / 3) * log((N / π) * sin(π * ℓ / N)) + c₁

    # The point of the layer: a bag built the ordinary way reaches the relations,
    # with ℓ, L and the cut count derived rather than passed.
    b = bag(
        entanglement_entropy(Region(1:32...)) => ring(32),
        entanglement_entropy(Region(5:20...)) => ring(16),
        CentralCharge => c,
    )
    rows = finite_size_entropy_report(b, PBC(N); c₁=c₁, atol=1e-12)
    @test length(rows) == 2
    @test all(r -> r.pass, rows)
    @test all(r -> r.relation isa CFTEntanglementPBC, rows)

    # A wrong entropy has to fail, or passing says nothing.
    bad = bag(entanglement_entropy(Region(1:32...)) => ring(32) + 0.1, CentralCharge => c)
    @test !only(finite_size_entropy_report(bad, PBC(N); c₁=c₁, atol=1e-12)).pass

    # Only the geometry each equation was derived for is matched, and the rest is
    # skipped rather than answered: four cuts on a ring, and no adjacency at all.
    mixed = bag(
        entanglement_entropy(Region(1, 2, 10, 11)) => 1.0,
        entanglement_entropy(Region("a", "b")) => 1.0,
        CentralCharge => c,
    )
    @test isempty(finite_size_entropy_report(mixed, PBC(N); c₁=c₁, atol=1e-12))

    # The same sixteen sites match at an open end and not in the bulk, which is
    # the distinction the cut count exists for.
    at_end = bag(entanglement_entropy(Region(1:16...)) => 0.0, CentralCharge => c)
    in_bulk = bag(entanglement_entropy(Region(20:35...)) => 0.0, CentralCharge => c)
    @test length(finite_size_entropy_report(at_end, OBC(N); c₁=c₁, atol=1e-12)) == 1
    @test only(finite_size_entropy_report(at_end, OBC(N); c₁=c₁, atol=1e-12)).relation isa
        CFTEntanglementOBC
    @test isempty(finite_size_entropy_report(in_bulk, OBC(N); c₁=c₁, atol=1e-12))

    # An infinite chain has no L to supply and gets Eq. (4).
    binf = bag(
        entanglement_entropy(Region(1:8...)) => (c / 3) * log(8) + c₁, CentralCharge => c
    )
    inf_rows = finite_size_entropy_report(binf, Infinite(); c₁=c₁, atol=1e-12)
    @test only(inf_rows).relation isa CFTEntanglementInfinite
    @test only(inf_rows).pass

    # Without a central charge there is nothing to check against, and no rows.
    @test isempty(
        finite_size_entropy_report(
            bag(entanglement_entropy(Region(1:32...)) => 1.0), PBC(N); c₁=c₁
        ),
    )

    # An empty report must mean "nothing matched", so a boundary condition with no
    # form is refused rather than quietly returning one.
    @test_throws "no finite-size form registered" finite_size_entropy_report(
        b, _NoFormBC(); c₁=c₁, atol=1e-12
    )
end

@testset "the random ring is reached the same way, with f sampled at ℓ/L" begin
    c̃, c₁′, N = log(2) / 2, 0.31, 1024
    conf(v) = sin(π * v) / π
    S̄(ℓ) = (c̃ / 3) * log(N * conf(ℓ / N)) + c₁′

    b = bag(entanglement_entropy(Region(1:300...)) => S̄(300), EffectiveCentralCharge => c̃)
    # `f` is a callable here, so the sweep samples it at the point the geometry
    # picks rather than trusting a number the caller computed.
    rows = finite_size_entropy_report(b, PBC(N); c₁=0.0, f=conf, c₁′=c₁′, atol=1e-12)
    @test only(rows).relation isa InfiniteRandomnessEntanglementPBC
    @test only(rows).pass

    # Without the scaling function the random form cannot be instantiated at all.
    @test isempty(finite_size_entropy_report(b, PBC(N); c₁=0.0, c₁′=c₁′, atol=1e-12))

    # A bag carrying both charges is checked against both laws, which is how the
    # clean and the random readings of one measurement are told apart.
    both = bag(
        entanglement_entropy(Region(1:300...)) => S̄(300),
        CentralCharge => 0.5,
        EffectiveCentralCharge => c̃,
    )
    kinds = [
        nameof(typeof(r.relation)) for
        r in finite_size_entropy_report(both, PBC(N); c₁=0.0, f=conf, c₁′=c₁′, atol=1e-12)
    ]
    @test Set(kinds) == Set([:CFTEntanglementPBC, :InfiniteRandomnessEntanglementPBC])
end

@testset "the scaling function is checked against what Eq. (24) says it is" begin
    # `f` reaches the relation as a number, so the relation can only ask that it be
    # positive. The callable is visible here, and Eq. (24) states it through an
    # expansion, `f(v) = Σₖ Aₖ sin((2k-1)πv)` under `Σₖ Aₖ(2k-1)π = 1`, from which
    # two properties follow without knowing a single coefficient.
    c̃, c₁′, N = log(2) / 2, 0.31, 1024
    conf(v) = sin(π * v) / π
    S̄(ℓ) = (c̃ / 3) * log(N * conf(ℓ / N)) + c₁′
    b = bag(entanglement_entropy(Region(1:300...)) => S̄(300), EffectiveCentralCharge => c̃)
    call(g) = finite_size_entropy_report(b, PBC(N); c₁=0.0, f=g, c₁′=c₁′, atol=1e-12)

    # The point of Eq. (24) is that the random ring is NOT the conformal one, and the
    # difference is the higher harmonics. So the check has to admit them: a two-term
    # f, normalised the same way, is exactly as legitimate and must pass. A check
    # that only accepted `sin(πv)/π` would refuse the physics it exists for.
    two(v) = (1 / π - 0.06) * sin(π * v) + 0.02 * sin(3π * v)
    @test isapprox(two(1e-5) / 1e-5, 1; atol=1e-6)      # the normalisation, restated
    @test !isapprox(two(0.3), conf(0.3); rtol=1e-3)     # and it really is a different f
    @test length(call(two)) == 1                        # accepted, and the row is real
    @test !only(call(two)).pass                         # against an entropy built on `conf`

    # Every basis function is symmetric about v = 1/2, so an admissible f is. An
    # asymmetric one is not in the family and is refused rather than scored.
    @test_throws "not reflection symmetric" call(v -> sin(π * v) / π + 0.05 * v)

    # `Σₖ Aₖ(2k-1)π = 1` IS `f'(0) = 1`. Off by a constant, `ln[L f]` shifts by `ln`
    # of it, which is absorbed by `c₁′` rather than showing up as a failure: the
    # wrong f would be reported as a pass with a wrong constant.
    @test_throws "f(v)/v" call(v -> sin(π * v))         # the chord, missing the 1/π
    @test_throws "f(v)/v" call(v -> 2 * sin(π * v) / π)

    # It refuses only where it is read. An `f` no row consumes changes nothing, and
    # refusing it would be refusing an argument for being present.
    open_chain = bag(
        entanglement_entropy(Region(1:300...)) => 1.0, EffectiveCentralCharge => c̃
    )
    @test isempty(
        finite_size_entropy_report(
            open_chain, OBC(N); c₁=0.0, f=v -> -1.0, c₁′=c₁′, atol=1e-12
        ),
    )
    clean = bag(entanglement_entropy(Region(1:300...)) => 1.0, CentralCharge => 0.5)
    @test length(
        finite_size_entropy_report(clean, PBC(N); c₁=0.0, f=v -> -1.0, atol=1e-12)
    ) == 1
end

@testset "an infinite chain gives the slope relations their derivative exactly" begin
    c, c₁ = 0.5, 0.4785
    S(ℓ) = (c / 3) * log(ℓ) + c₁
    b = bag(
        entanglement_entropy(Region(1:8...)) => S(8),
        entanglement_entropy(Region(1:32...)) => S(32),
        entanglement_entropy(Region(1:64...)) => S(64),
        CentralCharge => c,
    )
    rows = finite_size_entropy_report(b, Infinite(); c₁=c₁, atol=1e-12)
    slopes = filter(r -> r.relation isa CFTEntanglementSlope, rows)

    # `S` is affine in `ln ℓ` here, so the secant is the derivative and the rows
    # are exact rather than fitted: two consecutive pairs from three regions.
    @test length(slopes) == 2
    @test all(r -> r.pass, rows)
    @test all(r -> length(r.regions) == 2, slopes)

    # A finite chain's abscissa is the chord, so no slope row is taken there; the
    # closed form covers it and reporting an asymptotic slope as exact would not.
    bring = bag(
        entanglement_entropy(Region(1:8...)) => 0.0,
        entanglement_entropy(Region(1:16...)) => 0.0,
        CentralCharge => c,
    )
    @test !any(
        r -> r.relation isa CFTEntanglementSlope,
        finite_size_entropy_report(bring, PBC(64); c₁=c₁, atol=1e-12),
    )

    # The same measurement read at the other fixed point has to fail, or the two
    # readings are not being told apart.
    birfp = bag(
        entanglement_entropy(Region(1:8...)) => S(8),
        entanglement_entropy(Region(1:32...)) => S(32),
        EffectiveCentralCharge => log(2) / 2,
    )
    r = only(finite_size_entropy_report(birfp, Infinite(); c₁=0.0, atol=1e-12))
    @test r.relation isa InfiniteRandomnessEntanglementSlope
    @test !r.pass

    # And an actual random chain passes it, so the failure above is the reading
    # and not the relation.
    c̃ = log(2) / 2
    S̃(ℓ) = (c̃ / 3) * log(ℓ) + 0.31
    bok = bag(
        entanglement_entropy(Region(1:8...)) => S̃(8),
        entanglement_entropy(Region(1:32...)) => S̃(32),
        EffectiveCentralCharge => c̃,
    )
    @test only(finite_size_entropy_report(bok, Infinite(); c₁=0.0, atol=1e-12)).pass
end

@testset "a correlation length in the bag reaches the saturated form" begin
    c, ξ, N = 0.5, 12.0, 128
    b = bag(
        entanglement_entropy(Region(1:40...)) => 2 * (c / 6) * log(ξ),
        CentralCharge => c,
        CorrelationLength => ξ,
    )
    kinds = Dict(
        nameof(typeof(r.relation)) => r.pass for
        r in finite_size_entropy_report(b, PBC(N); c₁=0.4785, atol=1e-12)
    )

    # Off criticality the entropy sits on ξ, so the critical form fails and the
    # saturated one passes: the report says which regime the data is in.
    @test kinds[:OffCriticalEntanglementSaturation]
    @test !kinds[:CFTEntanglementPBC]

    # Only regions larger than ξ are matched, since below it there is no plateau
    # to be on.
    small = bag(
        entanglement_entropy(Region(1:8...)) => 0.0,
        CentralCharge => c,
        CorrelationLength => ξ,
    )
    @test !any(
        r -> r.relation isa OffCriticalEntanglementSaturation,
        finite_size_entropy_report(small, PBC(N); c₁=0.4785, atol=1e-12),
    )
end

@testset "the sweep refuses what it cannot pick, and what cannot be saturated" begin
    c, c₁, ξ = 0.5, 0.4785, 12.0

    # Two regions of one length carry one abscissa, and which a secant used would be
    # decided by Dict iteration order. That is a refusal, not a silent choice.
    tie = bag(
        entanglement_entropy(Region(1:4...)) => 0.5,
        entanglement_entropy(Region(101:104...)) => 0.9,
        entanglement_entropy(Region(1:16...)) => 1.1,
        CentralCharge => c,
    )
    @test_throws "no abscissa" finite_size_entropy_report(tie, Infinite(); c₁=c₁)

    # The saturated form needs both sides of the cut to be bulk, since the source
    # writes it for an infinite chain. A ring minus one site clears `ℓ > ξ` by a wide
    # margin, and its complement is a single site, which cannot hold the entropy the
    # row would otherwise pass.
    N = 2000
    almost_all = bag(
        entanglement_entropy(Region(setdiff(1:N, [1000])...)) => 2 * (c / 6) * log(ξ),
        CentralCharge => c,
        CorrelationLength => ξ,
    )
    @test !any(
        r -> r.relation isa OffCriticalEntanglementSaturation,
        finite_size_entropy_report(almost_all, PBC(N); c₁=c₁, atol=1e-12),
    )

    # A block with bulk on both sides is still matched, so the guard is not blanket.
    bulk = bag(
        entanglement_entropy(Region(1:200...)) => 2 * (c / 6) * log(ξ),
        CentralCharge => c,
        CorrelationLength => ξ,
    )
    @test any(
        r -> r.relation isa OffCriticalEntanglementSaturation && r.pass,
        finite_size_entropy_report(bulk, PBC(N); c₁=c₁, atol=1e-12),
    )
end
