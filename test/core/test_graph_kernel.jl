# The abstract typed-graph kernel (graph/kernel.jl) — the parent every concrete
# knowledge graph (quantity, derivation, model) instantiates.  Tested on a small
# synthetic graph over Symbol nodes, independent of any physics, so the generic
# traversal/reachability/path/export are verified in isolation.

using AbstractQAtlas
using AbstractQAtlas:
    TypedEdge,
    KnowledgeGraph,
    graph_edges,
    graph_nodes,
    graph_neighbors,
    graph_reachable,
    graph_shortest_path,
    graph_jsonl

@testset "TypedEdge fields and defaults" begin
    e = TypedEdge(:k, :a, :b, "lbl")
    @test e.kind === :k
    @test e.from === :a
    @test e.to === :b
    @test e.detail == "lbl"
    @test e.directed                       # default directed
    @test !TypedEdge(:k, :a, :b, "lbl", false).directed
    @test e isa TypedEdge{Symbol}
end

@testset "KnowledgeGraph is an edge collection" begin
    es = [TypedEdge(:e, :a, :b, "1"), TypedEdge(:e, :b, :c, "2")]
    g = KnowledgeGraph(es)
    @test length(g) == 2
    @test eltype(typeof(g)) === TypedEdge{Symbol}
    @test !isempty(g)
    @test collect(g) == es                 # iterates its edges
    @test graph_edges(g) === es
    @test isempty(KnowledgeGraph(TypedEdge{Symbol}[]))
    @test graph_nodes(g) == Set([:a, :b, :c])
end

@testset "graph_neighbors: all incident edges, direction ignored" begin
    g = KnowledgeGraph([
        TypedEdge(:e, :a, :b, "1"), TypedEdge(:e, :b, :c, "2"), TypedEdge(:e, :d, :b, "3")
    ])
    nb = graph_neighbors(g, :b)
    @test length(nb) == 3                   # b is from OR to in all three
    @test Set(e.detail for e in nb) == Set(["1", "2", "3"])
    @test isempty(graph_neighbors(g, :z))
end

@testset "graph_reachable respects per-edge direction" begin
    # directed chain a → b → c, plus a symmetric side link c — d
    g = KnowledgeGraph([
        TypedEdge(:e, :a, :b, "", true),
        TypedEdge(:e, :b, :c, "", true),
        TypedEdge(:e, :c, :d, "", false),
    ])
    @test graph_reachable(g, :a) == Set([:a, :b, :c, :d])   # forward + symmetric hop
    @test graph_reachable(g, :c) == Set([:c, :d])           # cannot go back up a directed edge
    @test graph_reachable(g, :d) == Set([:d, :c])           # symmetric edge is bidirectional
    @test graph_reachable(g, :z) == Set([:z])               # isolated node reaches only itself
end

@testset "graph_shortest_path: direction, empties, unreachable" begin
    g = KnowledgeGraph([
        TypedEdge(:step, :a, :b, "ab", true),
        TypedEdge(:step, :b, :c, "bc", true),
        TypedEdge(:link, :x, :y, "xy", false),
    ])
    # directed path follows the arrows
    p = graph_shortest_path(g, :a, :c)
    @test p !== nothing
    @test [e.detail for e in p] == ["ab", "bc"]
    # backward along directed edges is impossible
    @test graph_shortest_path(g, :c, :a) === nothing
    # symmetric edge traverses both ways
    @test length(graph_shortest_path(g, :y, :x)) == 1
    # same node ⇒ empty path
    @test graph_shortest_path(g, :a, :a) == TypedEdge{Symbol}[]
    # different components ⇒ nothing
    @test graph_shortest_path(g, :a, :x) === nothing
end

@testset "graph_jsonl export shape" begin
    g = KnowledgeGraph([TypedEdge(:k, :a, :b, "lbl", true)])
    io = IOBuffer()
    graph_jsonl(io, g)
    lines = split(strip(String(take!(io))), '\n')
    @test occursin("\"nodes\":2", lines[1])
    @test occursin("\"edges\":1", lines[1])
    @test occursin("\"kind\":\"k\"", lines[2])
    @test occursin("\"from\":\"a\"", lines[2])
    @test occursin("\"to\":\"b\"", lines[2])
    @test occursin("\"detail\":\"lbl\"", lines[2])
    @test occursin("\"directed\":true", lines[2])
    # custom node label
    io2 = IOBuffer()
    graph_jsonl(io2, g; nodelabel=x -> uppercase(string(x)))
    @test occursin("\"from\":\"A\"", String(take!(io2)))
end
