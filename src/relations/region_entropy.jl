# relations/region_entropy.jl: auto-discovery over the REGIONS present in a bag
# (design §5/§8b, Phase-2): the entropy inequalities, and the finite-size forms
# that read a critical chain's central charge off the same entries.
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

- **Maximum entropy**, for every region, when `local_dim` is given:
  `S(A) ≤ |A| · ln(local_dim)`.  This one is opt-in because a [`Region`](@ref) is
  a set of site labels with no Hilbert space attached, so the sweep cannot know
  `d`; omit it and an impossible entropy (5 nats on one qubit) produces no row.

  One uniform `d` for every site.  On a mixed lattice the bound is then only as
  tight as the value passed, and too generous a `d` MASKS a real violation —
  `S = 1.75` on two qubits fails at `local_dim = 2` and passes at `3`.  Verify
  uniformity before relying on a pass; a per-site mapping is not supported yet.

A negative (conditional) mutual information — a broken MPS/ED entanglement
calculation — is caught for whichever regions expose it.

Complementarity needs no relation of its own: for a bag in which `S(A∪B) = 0`,
Araki–Lieb already reads `0 ≥ |S(A) − S(B)|`, so a pure global state with
`S(A) ≠ S(B)` fails on the row that is already there.

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
function region_report(b::Bag; atol=0, local_dim::Union{Nothing,Int}=nothing)
    # Checked here rather than in the per-family helper: a bag with no entropy
    # family never reaches that helper, and an invalid `local_dim` would then be
    # indistinguishable from "no data yet".
    local_dim === nothing ||
        local_dim > 1 ||
        throw(ArgumentError("local_dim must be > 1; got $local_dim"))
    out = RegionReportRow[]
    for Q in _region_entropy_families(b)
        _region_report_family!(out, _region_entropies(b, Q); atol=atol, local_dim=local_dim)
    end
    return out
end
export region_report

# one family's sweep.  Split out of `region_report` so that adding a second
# entropy family cannot accidentally let regions from one family match unions
# from the other: the matcher only ever sees a single family's Dict.
function _region_report_family!(
    out::Vector{RegionReportRow},
    ents::AbstractDict;
    atol=0,
    local_dim::Union{Nothing,Int}=nothing,
)
    regions = sort!(collect(keys(ents)); by=r -> (length(r.sites), string(r)))
    # Maximum entropy `S(A) ≤ |A| ln d` — opt-in; see the docstring above for why.
    if local_dim !== nothing
        meb = MaxEntropyBound()
        for A in regions
            s = residual(meb; S=ents[A], log_d=length(A) * log(local_dim))
            push!(out, RegionReportRow(meb, (A,), s, _passes(meb, s, atol)))
        end
    end
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

`true` answers over the rows that were DISCOVERED, which for the maximum-entropy
family means the ones `local_dim` enabled.  Omitting it does not make that bound
pass — it makes it absent, and a bag that violates it can still answer `true` on
the strength of the inequalities that were checked.  Pass `local_dim` to include
it.
"""
function region_check_all(b::Bag; atol=0, local_dim::Union{Nothing,Int}=nothing)
    # reuse the shared "≥1 match, all pass" rule (interface.jl) so it can't drift
    return _all_passed(region_report(b; atol=atol, local_dim=local_dim))
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

# ─── Finite-size forms over the same region-keyed entries ──────────────────
#
# The inequalities above hold for any regions; these hold for one region at a
# time and only where the geometry matches the equation, so the matcher is a
# sweep with an admission test rather than a combination search.

"""
    RegionFiniteSizeRow

One row of a [`finite_size_entropy_report`](@ref): the `relation` it matched, the
`regions` it was auto-instantiated on (one for a closed form, the pair a slope
was taken across), its [`residual`](@ref), and `pass`.

Separate from [`RegionReportRow`](@ref) because the residual means something
else.  There it is a slack, satisfied at `≥ 0`; these are equalities, satisfied
only near `0`, and reading one as the other would call every negative residual a
violation and every large positive one a success.
"""
struct RegionFiniteSizeRow
    relation::AbstractRelation
    regions::Tuple{Vararg{Region}}
    residual::Number
    pass::Bool
end
export RegionFiniteSizeRow

# Sites `entanglement_cuts` can count on: nonempty, integer-labelled and not Bool,
# which is an Integer but not a lattice index. Shared, so the two sweeps below
# cannot admit a region the other skips.
function _countable_region(A::Region)
    return !isempty(A) && eltype(A.sites) <: Integer && eltype(A.sites) !== Bool
end

function _finite_size_row!(out, rel, regions::Tuple, vars, atol)
    r = residual(rel; vars...)
    push!(out, RegionFiniteSizeRow(rel, regions, r, _passes(rel, r, atol)))
    return out
end
function _finite_size_row!(out, rel, A::Region, vars, atol)
    return _finite_size_row!(out, rel, (A,), vars, atol)
end

# Two consequences of Eq. (24)'s `f(v) = Σₖ Aₖ sin((2k-1)πv)` under
# `Σₖ Aₖ(2k-1)π = 1`, neither needing a coefficient: the basis is symmetric about
# `v = 1/2`, and the normalisation is `f'(0) = 1`.  The FORM is not checked, the
# higher harmonics being exactly what a disorder average would be needed to see.
function _check_scaling_function(f)
    for v in (0.1, 0.25, 0.4)
        a, m = f(v), f(1 - v)
        isapprox(a, m; rtol=1e-8, atol=1e-12) || error(
            "finite_size_entropy_report: f($v) = $a but f($(1 - v)) = $m, so the " *
            "scaling function is not symmetric about v = 1/2 and is not in Eq. (24)'s " *
            "family.",
        )
    end
    v = 1e-5
    slope = f(v) / v
    isapprox(slope, 1; atol=1e-6) || error(
        "finite_size_entropy_report: f(v)/v → $slope, not 1, so the scaling function " *
        "breaks Eq. (24)'s normalisation. A constant factor shifts ln[L f] into c₁′ " *
        "and passes; f is the dimensionless sin(πv)/π, not the chord.",
    )
    return nothing
end

"""
    finite_size_entropy_report(b::Bag, bc::BoundaryCondition; c₁, kwargs...)
        -> Vector{RegionFiniteSizeRow}

Check every region entropy in `b` against the finite-size form its geometry
admits, under the boundary condition `bc`.

`ℓ` is the region's own length, `L` is `bc.N`, and the cut count comes from
[`entanglement_cuts`](@ref), so the three arguments most easily got wrong are
read off the bag rather than passed.  A region is matched only where its cut
count is the one its equation was derived for: two on a ring or an infinite
chain, one at an open end.  Anything else is skipped, as a non-disjoint pair is
skipped by [`region_report`](@ref), which is what makes a mixed bag usable.

On an infinite chain two regions also give the slope relations their derivative
exactly, since `S` is affine in `ln ℓ` there; a finite chain's abscissa is the
chord, so the closed forms cover it instead. [`HalvedChainEntropyDifference`](@ref)
is not reachable from one bag at all, needing entropies from two chain lengths.

The central charge is read from the bag, as `CentralCharge` and, when a random
critical chain is being checked, `EffectiveCentralCharge`; the latter also needs
`f`, the scaling function, and its own constant `c₁′`.  The non-universal
constants are arguments because they are not quantities: `c₁`, and `ln_g` for
the boundary entropy an open chain carries.

A region whose sites are not integers is skipped, having no adjacency to count
cuts with. A region of integer sites lying off the chain `bc` declares is not
skipped but refused, since that is the wrong `bc` for this bag and every later
row would be wrong the same way.

```julia
b = bag(entanglement_entropy(Region(1:32...)) => 1.06, CentralCharge => 0.5)
finite_size_entropy_report(b, PBC(64); c₁=0.4785)
```
"""
function finite_size_entropy_report(
    b::Bag, bc::BoundaryCondition; c₁::Real, ln_g::Real=0, f=nothing, c₁′::Real=0, atol=1e-8
)
    # An empty report must mean "no region matched", so a boundary condition with no
    # form here is refused rather than producing one: it would otherwise be a
    # boundary condition taught to `entanglement_cuts` and not to this sweep, and
    # the two are indistinguishable from the outside.
    bc isa Union{Infinite,OBC,PBC} || error(
        "finite_size_entropy_report: no finite-size form registered for $(typeof(bc))."
    )
    out = RegionFiniteSizeRow[]
    ents = _region_entropies(b)
    isempty(ents) && return out
    c = get(b, VariableKey(CentralCharge), nothing)
    c̃ = get(b, VariableKey(EffectiveCentralCharge), nothing)
    ξ = get(b, VariableKey(CorrelationLength), nothing)
    (c === nothing && c̃ === nothing) && return out
    # Only where it is read; an unconsumed `f` changes nothing.
    c̃ !== nothing && f !== nothing && bc isa PBC && _check_scaling_function(f)
    for pair in sort!(collect(ents); by=p -> (length(p.first), repr(p.first)))
        A, S = pair.first, pair.second
        (isempty(A) || !(eltype(A.sites) <: Integer) || eltype(A.sites) === Bool) &&
            continue
        ℓ, n = length(A), entanglement_cuts(bc, A)
        if c !== nothing
            if bc isa Infinite && n == 2
                _finite_size_row!(out, CFTEntanglementInfinite(), A, (; S, c, ℓ, c₁), atol)
            elseif bc isa PBC && n == 2
                _finite_size_row!(
                    out, CFTEntanglementPBC(), A, (; S, c, L=bc.N, ℓ, c₁), atol
                )
            elseif bc isa OBC && n == 1
                _finite_size_row!(
                    out, CFTEntanglementOBC(), A, (; S, c, L=bc.N, ℓ, c₁, ln_g), atol
                )
            end
        end
        if c̃ !== nothing && f !== nothing && bc isa PBC && n == 2
            _finite_size_row!(
                out,
                InfiniteRandomnessEntanglementPBC(),
                A,
                (; S̄=S, c̃, L=float(bc.N), f=f(ℓ / bc.N), c₁′),
                atol,
            )
        end
        # Off criticality the entropy stops following ℓ and sits on ξ instead, so
        # this is matched wherever a correlation length is in the bag and the region
        # is the larger of the two. Both sides of the cut, since the source writes it
        # as `S∞`: each cut saturates independently only if the complement is bulk
        # too, and a ring minus one site otherwise reports a pass on an entropy that
        # exceeds what a one-site subsystem can hold. The source states it for ξ ≪ ℓ,
        # and a row near ξ ≈ ℓ is expected to fail rather than be cut off by a
        # threshold chosen here.
        if c !== nothing && ξ !== nothing && ℓ > ξ && (bc isa Infinite || bc.N - ℓ > ξ)
            _finite_size_row!(
                out, OffCriticalEntanglementSaturation(), A, (; S, c, ξ, ncuts=n), atol
            )
        end
    end
    append!(out, _slope_rows(ents, bc, c, c̃, atol))
    return out
end

# The slope relations take a supplied derivative, and on an infinite chain two
# regions give it exactly rather than by fitting: `S` is affine in `ln ℓ` there, so
# the secant through any two points IS the derivative. On a finite chain the
# abscissa is the chord and not `ln ℓ`, so the closed forms above cover that case
# and this one stays out of it rather than reporting an asymptotic slope as exact.
function _slope_rows(ents, bc::BoundaryCondition, c, c̃, atol)
    out = RegionFiniteSizeRow[]
    bc isa Infinite || return out
    (c === nothing && c̃ === nothing) && return out
    pts = [
        (length(A), A, S) for (A, S) in ents if !isempty(A) &&
            eltype(A.sites) <: Integer &&
            eltype(A.sites) !== Bool &&
            entanglement_cuts(bc, A) == 2
    ]
    length(pts) >= 2 || return out
    sort!(pts; by=first)
    # Two regions of one length carry one abscissa, and which of them a secant used
    # would be decided by the order a Dict happened to iterate in. The law says their
    # entropies agree, so a bag holding both is either redundant or inconsistent, and
    # either way this cannot pick.
    for ((ℓ1, A1, _), (ℓ2, A2, _)) in zip(pts, pts[2:end])
        ℓ1 == ℓ2 && error(
            "finite_size_entropy_report: $A1 and $A2 both have length $ℓ1, so a slope " *
            "through them has no abscissa. Keep one, or give them distinct lengths.",
        )
    end
    for ((ℓ1, A1, S1), (ℓ2, A2, S2)) in zip(pts, pts[2:end])
        slope = (S2 - S1) / (log(ℓ2) - log(ℓ1))
        c === nothing || _finite_size_row!(
            out, CFTEntanglementSlope(), (A1, A2), (; dS_dlogℓ=slope, c, ncuts=2), atol
        )
        c̃ === nothing || _finite_size_row!(
            out,
            InfiniteRandomnessEntanglementSlope(),
            (A1, A2),
            (; dS_dlogℓ=slope, c̃, ncuts=2),
            atol,
        )
    end
    return out
end
export finite_size_entropy_report
