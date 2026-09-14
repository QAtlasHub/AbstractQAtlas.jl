# The registry as a directed derivation graph: reachability, lazy path-finding,
# and the traced solver.  Tests check that the graph is EQUALITIES-only (no
# inequality saturation leaks in as a "derivation"), that a found route
# actually computes the target (exactly, in rationals), that the debug trace
# names the true indirect route, and that unreachable targets fail loudly.

using AbstractQAtlas
using AbstractQAtlas:
    derive,
    derivable,
    derivation_steps,
    DerivationStep,
    DerivationTrace,
    AbstractInequality,
    solve,
    FreeEnergyFromZ,
    Rushbrooke,
    Widom

@testset "derivation graph is equalities-only (no inequality saturation)" begin
    steps = derivation_steps()
    @test !isempty(steps)
    # not one edge comes from an inequality — its `solve` returns a BOUND, not
    # an equational value, and must never masquerade as a derivation
    @test all(s -> !(s.relation isa AbstractInequality), steps)
    # every edge's inputs are exactly the relation's other variables
    for s in steps
        vs = AbstractQAtlas.variables(s.relation)
        @test s.output in vs
        @test Set(s.inputs) == Set(v for v in vs if v !== s.output)
    end
end

@testset "single-step derive matches a direct solve, with an indirect trace" begin
    v = derive(:f; Z=2.0, β=1.0)
    @test v ≈ solve(FreeEnergyFromZ(), Val(:f); Z=2.0, β=1.0)
    @test v ≈ -log(2.0)                       # F = −β⁻¹ ln Z
    t = derive(:f; Z=2.0, β=1.0, debug=true)
    @test t isa DerivationTrace
    @test t.indirect
    @test t.value ≈ v
    @test length(t.steps) == 1
    @test t.steps[1].relation isa FreeEnergyFromZ
    @test t.steps[1].output === :f
end

@testset "a directly-supplied target is flagged direct, no steps" begin
    t = derive(:Z; Z=5.0, β=1.0, debug=true)
    @test !t.indirect
    @test isempty(t.steps)
    @test t.value == 5.0
    @test derive(:Z; Z=5.0) == 5.0            # non-debug returns the bare value
end

@testset "multi-step route is found and computes exactly (rationals)" begin
    # {α, β} ─Rushbrooke→ γ ─Widom→ δ : two exact rational hops
    @test derive(:δ; α=0 // 1, β=1 // 8) == 15 // 1
    t = derive(:δ; α=0 // 1, β=1 // 8, debug=true)
    @test t.value == 15 // 1                   # Rational in ⇒ Rational out through the chain
    @test length(t.steps) == 2
    @test t.steps[1].relation isa Rushbrooke && t.steps[1].output === :γ
    @test t.steps[2].relation isa Widom && t.steps[2].output === :δ
    # dependency order: γ is produced before it is consumed
    @test findfirst(s -> s.output === :γ, t.steps) < findfirst(s -> :γ in s.inputs, t.steps)
end

@testset "the route is pruned to only contributing steps" begin
    # give extra, irrelevant knowns; the trace to :γ must not include them
    t = derive(:γ; α=0 // 1, β=1 // 8, δ=15 // 1, ν=1 // 1, η=1 // 4, debug=true)
    @test t.value == 7 // 4
    @test all(s -> s.output === :γ, t.steps)   # exactly the γ-producing step(s)
    @test length(t.steps) == 1
end

@testset "reachability: derivable is honest and self-consistent" begin
    r = derivable(; Z=2.0, β=1.0)
    @test :Z in r && :β in r                   # knowns are included
    @test :f in r                              # F reachable
    # everything reachable is actually derivable to a real value
    for sym in r
        @test derive(sym; Z=2.0, β=1.0) !== nothing
    end
    # a target outside the reachable set is NOT claimed
    @test :χ ∉ r
end

@testset "derivation_graph is a directed KnowledgeGraph instance" begin
    using AbstractQAtlas: derivation_graph, KnowledgeGraph, TypedEdge, graph_reachable
    dg = derivation_graph()
    @test dg isa KnowledgeGraph{Symbol}
    @test !isempty(dg)
    @test all(e -> e.directed, dg)          # every derivation edge is directed
    # each edge is input →[relation] output of some derivation step
    outs = Set(s.output for s in derivation_steps())
    ins = Set(i for s in derivation_steps() for i in s.inputs)
    for e in dg
        @test e.to in outs
        @test e.from in ins
    end
    # structural reachability OVER-approximates computability: the simple-edge
    # projection drops the hyperedges' AND-semantics.  Concretely, FreeEnergy
    # F = U − TS needs U AND S AND β, but the projection gives an edge U → F,
    # so `:U` structurally "reaches" `:F` while it cannot honestly derive it —
    # the whole point of keeping `derivable`/`derive` as the safe evaluators.
    @test :F in graph_reachable(dg, :U)          # structural: U alone reaches F
    @test :F ∉ derivable(; U=1.5)                # honest: cannot compute F from U alone
end

@testset "an unreachable target fails loudly (never a silent bad value)" begin
    err = try
        derive(:σxy; α=0 // 1)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("not reachable", err.msg)
end

@testset "every route back to a held-out exponent agrees, over the whole registry" begin
    # The hand-written cross-checks elsewhere pick their routes. This picks none:
    # each variable is held out in turn and every relation that reaches it from the
    # rest is run, so a wrong identity is caught wherever it sits rather than only
    # where someone thought to look.
    ising2d = (; α=0 // 1, β=1 // 8, γ=7 // 4, δ=15 // 1, ν=1 // 1, η=1 // 4, d=2 // 1)
    rows = consistency_report(ising2d; domain=:scaling)
    @test length(rows) == 7
    @test all(r -> r.agree, rows)
    # Rationals in, so the agreement is exact and not a tolerance being generous.
    @test all(r -> r.spread == 0, rows)
    @test consistent(ising2d; domain=:scaling)

    # γ is the busiest node, reached three ways; α and ν two ways each. A report
    # that found one route everywhere would be checking nothing.
    byname = Dict(r.target => length(r.steps) for r in rows)
    @test byname[:γ] == 3
    @test byname[:α] == 2
    @test byname[:ν] == 2
    @test Set(nameof(typeof(s.relation)) for r in rows for s in r.steps) ==
        Set([:Rushbrooke, :Widom, :Fisher, :Josephson])

    # A second consistent point, so the pass is not a property of one exponent set.
    mf = (; α=0 // 1, β=1 // 2, γ=1 // 1, δ=3 // 1, ν=1 // 2, η=0 // 1, d=4 // 1)
    @test consistent(mf; domain=:scaling)

    # And it discriminates: one exponent moved by 5% is denied through every route
    # that touches it, which is six of the seven rather than only γ's own.
    bad = (; α=0.0, β=0.125, γ=1.75 * 1.05, δ=15.0, ν=1.0, η=0.25, d=2.0)
    @test count(r -> !r.agree, consistency_report(bad; domain=:scaling)) == 6
    @test !consistent(bad; domain=:scaling)
end

@testset "a shared variable name is not a shared quantity" begin
    # The registry is one namespace and `S` is a ring's block in one relation and an
    # open chain's end block in another. One number cannot be both, so an unscoped
    # report says they disagree, correctly: the caller has supplied a state that does
    # not exist. `domain` is the remedy, and this pins that it is needed rather than
    # leaving a future reader to discover it from a puzzling row.
    c, c₁, N, ℓ = 0.5, 0.4785, 4096.0, 1024.0
    ring = (c / 3) * log((N / π) * sin(π * ℓ / N)) + c₁
    mixed = (; S=ring, c=c, L=N, ℓ=ℓ, c₁=c₁, ncuts=2, ln_g=0.0)
    row = only(filter(r -> r.target === :S, consistency_report(mixed)))
    @test length(row.steps) >= 2
    @test !row.agree

    # Scoped to the one geometry the number belongs to, the same data is consistent.
    @test solve(CFTEntanglementPBC(), Val(:S); c=c, L=N, ℓ=ℓ, c₁=c₁) ≈ ring atol = 1e-12

    # Nothing reachable is reported as agreement, so `consistent` alone is not proof
    # that anything was checked.
    @test isempty(consistency_report((; α=0.11)))
    @test consistent((; α=0.11))
end

@testset "applicability is not inferred, and the report says where it stopped" begin
    # A classical exponent set agrees over `:scaling`. It does so because the
    # relations that do not apply at a classical point happened to need an input
    # that was absent, not because anything checked applicability.
    classical = (; α=0 // 1, β=1 // 8, γ=7 // 4, δ=15 // 1, ν=1 // 1, η=1 // 4, d=2 // 1)
    @test consistent(classical; domain=:scaling)

    # Supplying `z`, which a classical point may perfectly well have, lets
    # QuantumHyperscaling reach α from the same knowns, and four rows break. The
    # identities are both right; `2 - α = (d + z)ν` simply is not a law at this
    # point, and nothing in the registry says so.
    withz = (; classical..., z=1 // 1)
    broken = consistency_report(withz; domain=:scaling)
    @test !all(r -> r.agree, broken)
    α_row = only(filter(r -> r.target === :α, broken))
    @test :QuantumHyperscaling in Set(nameof(typeof(s.relation)) for s in α_row.steps)
    @test :Josephson in Set(nameof(typeof(s.relation)) for s in α_row.steps)

    # `exclude` is how a caller states the scope, and with the quantum law out the
    # same data is consistent again. That the answer turns on one exclusion is the
    # point: the tool checks algebra, and the physics of which laws hold together
    # is the caller's to declare.
    @test consistent(withz; domain=:scaling, exclude=(:QuantumHyperscaling,))

    # Excluding the wrong one does not rescue it, so the knob is not a way to make
    # any data pass.
    @test !consistent(withz; domain=:scaling, exclude=(:Fisher,))
end

@testset "the type-keyed cross-check is sound where the name-keyed one cannot be" begin
    # A VariableKey is a quantity type, so two relations meet at a node only when
    # they mean the same thing. The name-keyed graph cannot do this: `β` is the
    # order-parameter exponent in Rushbrooke and the inverse temperature in
    # DetailedBalance, one symbol produced by nineteen relations.
    produce_β = unique(
        nameof(typeof(st.relation)) for
        st in AbstractQAtlas.derivation_steps() if st.output === :β
    )
    @test length(produce_β) >= 15
    @test :Rushbrooke in produce_β
    @test :DetailedBalance in produce_β        # inverse temperature, not an exponent

    # Typed, those are separate nodes and cannot be compared with each other.
    @test VariableKey(InverseTemperature) != VariableKey(CentralCharge)
    ising_bag = bag(
        SpecificHeatExponent => 0 // 1,
        OrderParameterExponent => 1 // 8,
        SusceptibilityExponent => 7 // 4,
        CriticalIsothermExponent => 15 // 1,
        CorrelationLengthExponent => 1 // 1,
        AnomalousDimension => 1 // 4,
        SpatialDimension => 2 // 1,
    )
    rows = consistency_report(ising_bag)
    @test length(rows) == 7
    @test all(r -> r.agree, rows)
    @test all(r -> r.spread == 0, rows)
    @test Set(nameof(r.target.type) for r in rows) == Set([
        :SpecificHeatExponent,
        :OrderParameterExponent,
        :SusceptibilityExponent,
        :CriticalIsothermExponent,
        :CorrelationLengthExponent,
        :AnomalousDimension,
        :SpatialDimension,
    ])

    # The classical identities reach the typed graph only because their exponents
    # were given quantities of their own in this branch. Before that `Rushbrooke`
    # had no typed slot and was absent here, which is what left `:β` the exponent
    # and `:β` the inverse temperature sharing a node on the name-keyed side.
    @test variable_types(Rushbrooke()) ==
        (SpecificHeatExponent, OrderParameterExponent, SusceptibilityExponent)
    @test any(st -> st.relation isa Rushbrooke, AbstractQAtlas.typed_derivation_steps())
    @test any(r -> :Rushbrooke in Set(nameof(typeof(s.relation)) for s in r.steps), rows)

    # The name-keyed report still reaches them too, and agrees, so typing narrowed
    # nothing: the two now overlap on this family instead of being disjoint.
    ising2d = (; α=0 // 1, β=1 // 8, γ=7 // 4, δ=15 // 1, ν=1 // 1, η=1 // 4, d=2 // 1)
    @test consistent(ising2d; domain=:scaling)

    # Typing does not settle everything: `d` is one SpatialDimension by design, and
    # the classical image's dimension and the chain's own are the same node, so the
    # reading a caller means is still theirs to state.
    @test SpatialDimension in variable_types(Josephson())
    @test SpatialDimension in variable_types(QuantumHyperscaling())
end

@testset "alternatives are grouped, so a law that does not apply is not a contradiction" begin
    # Iglói and Monthus state one menu of observables per fixed-point type, so the
    # conventional and the activated susceptibility are two readings of one quantity
    # and never both true. `law_family` records which compete.
    @test law_family(Josephson()) === law_family(QuantumHyperscaling())
    @test law_family(ConventionalFieldSusceptibility()) ===
        law_family(ActivatedSusceptibility()) ===
        law_family(GriffithsSusceptibility())
    @test law_family(ConventionalFiniteSizeEnergy()) ===
        law_family(ActivatedFiniteSizeScaling())
    # A law with no alternative is alone, so nothing that had no family gained one.
    @test law_family(Rushbrooke()) === :Rushbrooke
    @test law_family(Widom()) !== law_family(Fisher())

    classical = (; α=0 // 1, β=1 // 8, γ=7 // 4, δ=15 // 1, ν=1 // 1, η=1 // 4, d=2 // 1)
    withz = (; classical..., z=1 // 1)
    bad = filter(r -> !r.agree, consistency_report(withz; domain=:scaling))

    # Before grouping this broke four rows. Three of them had an applicable
    # alternative in the family, Josephson beside QuantumHyperscaling, and a family
    # holds when one member does.
    @test length(bad) == 1
    @test only(bad).target === :z

    # The one that remains is the boundary: `z` is reached only by the law that does
    # not apply, so its family has a single member and there is nothing to stand in.
    # No grouping fixes that; `exclude` is the caller saying so.
    @test Set(nameof(typeof(s.relation)) for s in only(bad).steps) ==
        Set([:QuantumHyperscaling])
    @test consistent(withz; domain=:scaling, exclude=(:QuantumHyperscaling,))

    # Grouping must not soften the check. A wrong exponent is denied as widely as
    # before, and a family whose members all miss is still caught.
    wrongγ = (; α=0.0, β=0.125, γ=1.75 * 1.05, δ=15.0, ν=1.0, η=0.25, d=2.0)
    @test count(r -> !r.agree, consistency_report(wrongγ; domain=:scaling)) == 6
    wrongα = (; α=0.9, β=0.125, γ=1.75, δ=15.0, ν=1.0, η=0.25, d=2.0, z=1.0)
    @test !consistent(wrongα; domain=:scaling)
end

@testset "no untyped symbol carries two meanings across the registry" begin
    # The name-keyed graph makes one node per SYMBOL, so a letter that means two
    # quantities is a node that compares them. `β` was the exemplar, the exponent
    # against the inverse temperature; typing the scaling side split it. This is the
    # standing sweep that keeps the next one from arriving unnoticed: a symbol
    # appearing across three or more domains with no relation typing it anywhere.
    bysym = Dict{Symbol,Set{Symbol}}()
    typedsym = Dict{Symbol,Set{Any}}()
    for r in all_relations()
        d = AbstractQAtlas.domain(r)
        slots = Dict(variable_slots(r))
        for v in variables(r)
            push!(get!(bysym, v, Set{Symbol}()), d)
            k = get(slots, v, nothing)
            k === nothing || push!(get!(typedsym, v, Set{Any}()), k)
        end
    end
    bare = Set(
        s for (s, ds) in bysym if length(ds) >= 3 && isempty(get(typedsym, s, Set()))
    )

    # `ω` is a frequency in every one of its domains, so one node is the right
    # answer there. An allow-list rather than a threshold, so a new letter is
    # refused until someone says which quantity it is.
    @test bare == Set([:ω])

    # `ζ` was the one this sweep found: an eigenvalue of a Gaussian correlation
    # matrix in `EntanglementSpectrumCorrelation` and the large-spin exponent in
    # `LargeSpinMoment`, two quantities with nothing to do with each other.
    @test LargeSpinExponent in variable_types(LargeSpinMoment())
    @test CorrelationMatrixEigenvalue in variable_types(EntanglementSpectrumCorrelation())
    @test LargeSpinExponent !== CorrelationMatrixEigenvalue

    # And the exemplar stays split.
    @test OrderParameterExponent in variable_types(Rushbrooke())
    @test !(InverseTemperature in variable_types(Rushbrooke()))
end
