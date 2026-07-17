# relations/region_entropy.jl — auto-discovery of the entanglement-entropy
# inequalities over the REGIONS present in a bag (design §5/§8b, Phase-2).
#
# The entropy inequalities hold for ANY (disjoint) regions.  Keyed on a Region
# support (`entanglement_entropy(A)`), they become auto-discoverable: `region_report`
# scans a bag of region-entropies and checks, on every matching region combination,
# subadditivity + Araki–Lieb (disjoint PAIRS) and strong subadditivity (pairwise-
# disjoint TRIPLES) — no hand-labeled A/B/AB/ABC.  The relations' scalar kernels
# (Subadditivity, ArakiLieb, StrongSubadditivity) are reused verbatim; this is the
# region-matching layer over them.  Multipartite (KP/LW TEE) + §8b index-unification
# follow on this same matcher.

"""
    RegionReportRow

One row of a [`region_report`](@ref): the `relation` (an entropy inequality), the
pairwise-disjoint `regions` it was auto-instantiated on (`(A, B)` for the bipartite
inequalities, `(A, B, C)` for strong subadditivity), the `slack` (its
[`residual`](@ref); `≥ 0` ⇔ satisfied), and `pass`.
"""
struct RegionReportRow
    relation::AbstractRelation
    regions::Tuple{Vararg{Region}}
    slack::Number
    pass::Bool
end
export RegionReportRow

# the von Neumann region-entropies present in a bag: Region → S(Region)
function _region_entropies(b::Bag)
    return Dict(
        k.support.region => v for
        (k, v) in b if k.type === VonNeumannEntropy && k.support isa RegionSupport
    )
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

A negative (conditional) mutual information — a broken MPS/ED entanglement
calculation — is caught for whichever regions expose it.

```julia
b = bag(entanglement_entropy(1) => 0.7, entanglement_entropy(2) => 0.7,
        entanglement_entropy(1, 2) => 1.0)      # S(A), S(B), S(A∪B)
all(row -> row.pass, region_report(b))          # true — S is subadditive here
```
"""
function region_report(b::Bag; atol=0)
    ents = _region_entropies(b)
    regions = collect(keys(ents))
    out = RegionReportRow[]
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
    # strong subadditivity over pairwise-disjoint triples (A, B, C): B is the shared
    # middle, {A, C} unordered (SSA is symmetric in A↔C).
    for bi in eachindex(regions)
        B = regions[bi]
        for i in eachindex(regions), k in (i + 1):lastindex(regions)
            (i == bi || k == bi) && continue
            A, C = regions[i], regions[k]
            (disjoint(A, B) && disjoint(B, C) && disjoint(A, C)) || continue
            AB, BC, ABC = A ∪ B, B ∪ C, A ∪ B ∪ C
            (haskey(ents, AB) && haskey(ents, BC) && haskey(ents, ABC)) || continue
            rel = StrongSubadditivity()
            s = residual(rel; S_AB=ents[AB], S_BC=ents[BC], S_ABC=ents[ABC], S_B=ents[B])
            push!(out, RegionReportRow(rel, (A, B, C), s, _passes(rel, s, atol)))
        end
    end
    return out
end
export region_report

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
    mutual_information(b::Bag, A::Region, B::Region) -> Number

The mutual information `I(A:B) = S(A) + S(B) − S(A∪B)`, computed from the region
entropies in the bag `b` (the [`Subadditivity`](@ref) slack; `≥ 0`).  Errors if any
of the three entropies is absent.

```julia
mutual_information(bag(entanglement_entropy(1) => 0.7, entanglement_entropy(2) => 0.7,
                       entanglement_entropy(1, 2) => 1.0), Region(1), Region(2))   # 0.4
```
"""
function mutual_information(b::Bag, A::Region, B::Region)
    ents = _region_entropies(b)
    for R in (A, B, A ∪ B)
        haskey(ents, R) || error("mutual_information: S($R) is not in the bag")
    end
    return ents[A] + ents[B] - ents[A ∪ B]
end
export mutual_information
