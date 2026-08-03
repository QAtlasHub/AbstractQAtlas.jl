# relations/region_entropy.jl — auto-discovery of the entanglement-entropy
# inequalities over the REGIONS present in a bag (design §5/§8b, Phase-2).
#
# The entropy inequalities hold for ANY (disjoint) regions.  Keyed on a Region
# support (`entanglement_entropy(A)`), they become auto-discoverable: `region_report`
# scans a bag of region-entropies and checks, on every matching region combination,
# subadditivity + Araki–Lieb (disjoint PAIRS) and strong subadditivity (pairwise-
# disjoint TRIPLES) — no hand-labeled A/B/AB/ABC.  The relations' scalar kernels
# (Subadditivity, ArakiLieb, StrongSubadditivity) are reused verbatim; this is the
# region-matching layer over them.  Kitaev–Preskill TEE auto-discovery (`region_tee_report`)
# rides the same matcher; Levin–Wen TEE + §8b index-unification follow.

"""
    RegionReportRow

One row of a [`region_report`](@ref): the `relation` (an entropy inequality), the
pairwise-disjoint `regions` it was auto-instantiated on (`(A, B)` for the bipartite
inequalities, `(A, B, C)` for the triple ones — strong subadditivity and weak
monotonicity), the `slack` (its [`residual`](@ref); `≥ 0` ⇔ satisfied), and `pass`.
"""
struct RegionReportRow
    relation::AbstractRelation
    regions::Tuple{Vararg{Region}}
    slack::Number
    pass::Bool
end
export RegionReportRow

"""
    obeys_entropy_inequalities(::Type) -> Bool

Whether a region-keyed quantity satisfies subadditivity, Araki–Lieb, strong
subadditivity and weak monotonicity, and may therefore be swept by
[`region_report`](@ref).

Opt-in, defaulting to `false`, because "is an entanglement measure" is not the
criterion: [`RenyiEntropy`](@ref) and [`TsallisEntropy`](@ref) live on the same
regions and are **not** strongly subadditive away from the von Neumann limit, so
a `<: AbstractEntanglementMeasure` test would auto-discover inequalities they are
not required to satisfy and report correct data as broken.  Declared `true` only
for [`VonNeumannEntropy`](@ref) and [`FermionicEntanglementEntropy`](@ref), each
of which is the von Neumann entropy of an honest reduced state.
"""
obeys_entropy_inequalities(::Type) = false
obeys_entropy_inequalities(::Type{VonNeumannEntropy}) = true
obeys_entropy_inequalities(::Type{FermionicEntanglementEntropy}) = true
export obeys_entropy_inequalities

# the region-entropies of ONE quantity present in a bag: Region → S(Region)
#
# `===` and not `<:`, and the quantity is an ARGUMENT rather than hard-coded.  A
# bag may legitimately carry `VonNeumannEntropy` and `FermionicEntanglementEntropy`
# on the same regions holding DIFFERENT numbers — they agree only on a single
# contiguous interval — so merging the two families would build inequalities out of
# one entropy's `S(A)` and the other's `S(A∪B)`.  That mixture is not a quantity at
# all, and it fails or passes for no stateable reason.
#
# (An earlier comment here justified `<:` by `VonNeumannEntropy` being a parametric
# family with a `{:quench}` member.  It is not: it is a plain struct, and the quench
# entropy was split off downstream as `QuenchEntanglementEntropy`.  With no subtypes
# the two spellings agreed, so nothing failed when the rationale expired.)
function _region_entropies(b::Bag, Q::Type=VonNeumannEntropy)
    return Dict(
        k.support.region => v for (k, v) in b if k.type === Q && k.support isa RegionSupport
    )
end

# every entropy family in the bag that has opted into the inequalities, in a
# deterministic order so report rows do not shuffle between runs
function _region_entropy_families(b::Bag)
    ts = unique(
        k.type for
        (k, _) in b if k.support isa RegionSupport && obeys_entropy_inequalities(k.type)
    )
    return sort!(collect(ts); by=string)
end

"""
    region_report(b::Bag; atol=0) -> Vector{RegionReportRow}

Auto-discover the entanglement-entropy inequalities over the REGIONS in a bag of
region-keyed entropies (`bag(entanglement_entropy(A) => s_A, …)`), with no A/B/AB
hand-labeling — the region twin of [`relation_report`](@ref):

- **Subadditivity** and **Araki–Lieb**, for every disjoint pair `(A, B)` whose
  `S(A)`, `S(B)`, `S(A∪B)` are all present: `I(A:B) = S(A)+S(B)−S(A∪B) ≥ 0` and
  `S(A∪B) ≥ |S(A)−S(B)|`.
- **Strong subadditivity**, for every pairwise-disjoint triple `(A, B, C)` whose
  `S(B)`, `S(A∪B)`, `S(B∪C)`, `S(A∪B∪C)` are present:
  `S(A∪B) + S(B∪C) ≥ S(A∪B∪C) + S(B)` (the conditional mutual information
  `I(A:C|B) ≥ 0`).
- **Weak monotonicity**, for every pairwise-disjoint triple `(A, B, C)` whose
  `S(A)`, `S(C)`, `S(A∪B)`, `S(B∪C)` are present — no full-system `S(A∪B∪C)`, so it is
  found strictly more often than strong subadditivity: `S(A∪B) + S(B∪C) ≥ S(A) + S(C)`.

A negative (conditional) mutual information — a broken MPS/ED entanglement
calculation — is caught for whichever regions expose it.

Every entropy family in the bag that declares
[`obeys_entropy_inequalities`](@ref) is swept **separately**: a bag holding both
`entanglement_entropy(A)` and [`fermionic_entanglement_entropy`](@ref)`(A)` on
the same regions yields both sets of rows, and no inequality is ever built from
one family's `S(A)` and another's `S(A∪B)`.

```julia
b = bag(entanglement_entropy(1) => 0.7, entanglement_entropy(2) => 0.7,
        entanglement_entropy(1, 2) => 1.0)      # S(A), S(B), S(A∪B)
all(row -> row.pass, region_report(b))          # true — S is subadditive here
```
"""
function region_report(b::Bag; atol=0)
    out = RegionReportRow[]
    for Q in _region_entropy_families(b)
        _region_report_family!(out, _region_entropies(b, Q); atol=atol)
    end
    return out
end
export region_report

# one family's sweep.  Split out of `region_report` so that adding a second
# entropy family cannot accidentally let regions from one family match unions
# from the other: the matcher only ever sees a single family's Dict.
function _region_report_family!(out::Vector{RegionReportRow}, ents::AbstractDict; atol=0)
    regions = sort!(collect(keys(ents)); by=r -> (length(r.sites), string(r)))
    for i in eachindex(regions), j in (i + 1):lastindex(regions)
        A, B = regions[i], regions[j]
        disjoint(A, B) || continue
        haskey(ents, A ∪ B) || continue
        S_A, S_B, S_AB = ents[A], ents[B], ents[A ∪ B]
        for rel in (Subadditivity(), ArakiLieb())
            s = residual(rel; S_A=S_A, S_B=S_B, S_AB=S_AB)
            push!(out, RegionReportRow(rel, (A, B), s, _passes(rel, s, atol)))
        end
    end
    # strong subadditivity + weak monotonicity over pairwise-disjoint triples (A, B, C):
    # B is the shared middle, {A, C} unordered (both are symmetric in A↔C). Weak
    # monotonicity S(A∪B)+S(B∪C) ≥ S(A)+S(C) needs no full-system S(A∪B∪C), so it is
    # discovered whenever the two pair-unions are present — strictly more often than SSA.
    for bi in eachindex(regions)
        B = regions[bi]
        for i in eachindex(regions), k in (i + 1):lastindex(regions)
            (i == bi || k == bi) && continue
            A, C = regions[i], regions[k]
            (disjoint(A, B) && disjoint(B, C) && disjoint(A, C)) || continue
            AB, BC = A ∪ B, B ∪ C
            (haskey(ents, AB) && haskey(ents, BC)) || continue
            wm = WeakMonotonicity()
            sw = residual(wm; S_AB=ents[AB], S_BC=ents[BC], S_A=ents[A], S_C=ents[C])
            push!(out, RegionReportRow(wm, (A, B, C), sw, _passes(wm, sw, atol)))
            # strong subadditivity additionally needs the full-system entropy S(A∪B∪C)
            ABC = A ∪ B ∪ C
            haskey(ents, ABC) || continue
            ssa = StrongSubadditivity()
            s = residual(ssa; S_AB=ents[AB], S_BC=ents[BC], S_ABC=ents[ABC], S_B=ents[B])
            push!(out, RegionReportRow(ssa, (A, B, C), s, _passes(ssa, s, atol)))
        end
    end
    return out
end

"""
    region_check_all(b::Bag; atol=0) -> Bool

`true` iff every entropy inequality (bipartite + strong subadditivity) auto-discovered
by [`region_report`](@ref) holds on the bag `b` — and at least one instance was found
(an empty match is `false`, never a silent green).
"""
function region_check_all(b::Bag; atol=0)
    # reuse the shared "≥1 match, all pass" rule (interface.jl) so it can't drift
    return _all_passed(region_report(b; atol=atol))
end
export region_check_all

"""
    mutual_information(b::Bag, A::Region, B::Region; quantity=VonNeumannEntropy) -> Number

The mutual information `I(A:B) = S(A) + S(B) − S(A∪B)`, computed from the region
entropies in the bag `b` (the [`Subadditivity`](@ref) slack; `≥ 0`).  Errors if any
of the three entropies is absent.

`quantity` selects the entropy family — pass
[`FermionicEntanglementEntropy`](@ref) to read the fermionic mutual information
out of a bag built with [`fermionic_entanglement_entropy`](@ref).  The two are
different numbers whenever a region is disconnected, and the difference does not
cancel here, so the family is named rather than inferred.

```julia
mutual_information(bag(entanglement_entropy(1) => 0.7, entanglement_entropy(2) => 0.7,
                       entanglement_entropy(1, 2) => 1.0), Region(1), Region(2))   # 0.4
```
"""
function mutual_information(b::Bag, A::Region, B::Region; quantity=VonNeumannEntropy)
    ents = _region_entropies(b, quantity)
    for R in (A, B, A ∪ B)
        haskey(ents, R) || error("mutual_information: S($R) is not in the bag")
    end
    return ents[A] + ents[B] - ents[A ∪ B]
end
export mutual_information

# fetch S(R) for each region, erroring by name (`what`) if any is absent
function _region_S(b::Bag, what::String, regions...)
    ents = _region_entropies(b)
    for R in regions
        haskey(ents, R) || error("$what: S($R) is not in the bag")
    end
    return (ents[R] for R in regions)
end

# The multipartite combinations below are the named invariants they claim to be only on
# a genuine tripartition (pairwise-disjoint A, B, C) — the same precondition `region_report`
# enforces before auto-discovering SSA.  With an overlapping/repeated region the unions
# collapse (e.g. C == A ⇒ A∪B∪C = A∪B) and the sum silently returns a physical-looking but
# meaningless number, so guard it rather than trust the caller.
function _require_tripartition(what::String, A::Region, B::Region, C::Region)
    (disjoint(A, B) && disjoint(B, C) && disjoint(A, C)) ||
        error("$what: A, B, C must be pairwise disjoint")
    return nothing
end

"""
    conditional_mutual_information(b::Bag, A::Region, B::Region, C::Region) -> Number

The conditional mutual information
`I(A:C|B) = S(A∪B) + S(B∪C) − S(A∪B∪C) − S(B)`, computed from the region entropies
in `b` for pairwise-disjoint `A, B, C` (the [`StrongSubadditivity`](@ref) /
[`MarkovEntropyDefinition`](@ref) slack; `≥ 0` by SSA).  Errors if the regions are not
a tripartition or if any of the four entropies is absent.
"""
function conditional_mutual_information(b::Bag, A::Region, B::Region, C::Region)
    _require_tripartition("conditional_mutual_information", A, B, C)
    S_AB, S_BC, S_ABC, S_B = _region_S(
        b, "conditional_mutual_information", A ∪ B, B ∪ C, A ∪ B ∪ C, B
    )
    return S_AB + S_BC - S_ABC - S_B
end
export conditional_mutual_information

"""
    tripartite_information(b::Bag, A::Region, B::Region, C::Region) -> Number

The tripartite (interaction) information
`I₃ = S(A)+S(B)+S(C) − S(A∪B)−S(A∪C)−S(B∪C) + S(A∪B∪C) = I(A:B) + I(A:C) − I(A:B∪C)`,
from the region entropies in `b` for pairwise-disjoint `A, B, C` — equal to
`−`[`topological_entanglement_entropy`](@ref) (the Kitaev–Preskill combination).  Errors
if the regions are not a tripartition or if any of the seven entropies is absent.
"""
function tripartite_information(b::Bag, A::Region, B::Region, C::Region)
    _require_tripartition("tripartite_information", A, B, C)
    S_A, S_B, S_C, S_AB, S_AC, S_BC, S_ABC = _region_S(
        b, "tripartite_information", A, B, C, A ∪ B, A ∪ C, B ∪ C, A ∪ B ∪ C
    )
    return S_A + S_B + S_C - S_AB - S_AC - S_BC + S_ABC
end
export tripartite_information

"""
    topological_entanglement_entropy(b::Bag, A::Region, B::Region, C::Region) -> Number

The Kitaev–Preskill topological entanglement entropy `γ = ln 𝒟` from a tripartition
(Kitaev & Preskill, [KitaevPreskill2006](@cite)),
`γ = −[S(A)+S(B)+S(C) − S(A∪B)−S(B∪C)−S(C∪A) + S(A∪B∪C)]` — the area-law-independent
constant isolated by the alternating tripartite sum ([`KitaevPreskillTEE`](@ref);
`γ > 0` ⇒ topological order).  Equals `−`[`tripartite_information`](@ref).
"""
function topological_entanglement_entropy(b::Bag, A::Region, B::Region, C::Region)
    return -tripartite_information(b, A, B, C)
end
export topological_entanglement_entropy

"""
    RegionTEERow

One row of a [`region_tee_report`](@ref): the pairwise-disjoint tripartition `regions`
`(A, B, C)` it was auto-instantiated on, the tripartite information
`tripartite_information` (`I₃`), and the Kitaev–Preskill topological entanglement entropy
`topological_entanglement_entropy` (`γ = −I₃`).
"""
struct RegionTEERow
    regions::NTuple{3,Region}
    tripartite_information::Number
    topological_entanglement_entropy::Number
end
export RegionTEERow

"""
    region_tee_report(b::Bag) -> Vector{RegionTEERow}

Auto-discover the tripartite information `I₃` and the Kitaev–Preskill topological
entanglement entropy `γ = −I₃` over the REGIONS in a bag of region-keyed entropies — the
multipartite twin of [`region_report`](@ref) (which handles the entropy *inequalities*).
One row is emitted per pairwise-disjoint triple `{A, B, C}` whose seven sub-entropies
`S(A)`, `S(B)`, `S(C)`, `S(A∪B)`, `S(A∪C)`, `S(B∪C)`, `S(A∪B∪C)` are all present; `I₃` is
symmetric in `A, B, C`, so each unordered triple gives exactly one row.

`γ` is the [`KitaevPreskillTEE`](@ref) constant `ln 𝒟` — *provided the regions form a KP
tripartition* (three sectors meeting so the boundary-law terms cancel). The set layer
carries no geometry, so this reports the alternating sum for any admissible triple; whether
it isolates the topological constant is the caller's (geometry-dependent) responsibility.

```julia
γ = log(2)
b = bag(entanglement_entropy(1) => 1.0, entanglement_entropy(2) => 1.0,
        entanglement_entropy(3) => 1.0, entanglement_entropy(1, 2) => 1.5,
        entanglement_entropy(1, 3) => 1.5, entanglement_entropy(2, 3) => 1.5,
        entanglement_entropy(1, 2, 3) => 1.5 - γ)   # area terms cancel, leaving −γ
only(region_tee_report(b)).topological_entanglement_entropy ≈ γ   # ln 2 (toric code)
```
"""
function region_tee_report(b::Bag)
    ents = _region_entropies(b)
    regions = collect(keys(ents))
    out = RegionTEERow[]
    for i in eachindex(regions),
        j in (i + 1):lastindex(regions),
        k in (j + 1):lastindex(regions)

        A, B, C = regions[i], regions[j], regions[k]
        (disjoint(A, B) && disjoint(A, C) && disjoint(B, C)) || continue
        (
            haskey(ents, A ∪ B) &&
            haskey(ents, A ∪ C) &&
            haskey(ents, B ∪ C) &&
            haskey(ents, A ∪ B ∪ C)
        ) || continue
        # reuse the verified I₃ combination (single source of truth); the disjoint + haskey
        # gates above guarantee the helper's own guards pass, so it never errors here.
        I3 = tripartite_information(b, A, B, C)
        push!(out, RegionTEERow((A, B, C), I3, -I3))
    end
    return out
end
export region_tee_report
