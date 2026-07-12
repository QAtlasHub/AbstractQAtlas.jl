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
