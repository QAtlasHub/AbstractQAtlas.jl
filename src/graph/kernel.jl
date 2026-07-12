# graph/kernel.jl — the abstract typed-graph parent.
#
# One abstraction under which every knowledge graph in the ecosystem is an
# instance: the quantity-relationship graph (`structure/graph.jl`), the
# directed derivation graph (`relations/derivation.jl`), and — once QAtlas is
# refactored onto AbstractQAtlas — its model graph (realizes / reduces / dual /
# limits).  They share one shape: a bag of TYPED EDGES over nodes of some type
# `N` (a quantity `Type`, a relation variable `Symbol`, a model `Type`).  The
# generic traversal / reachability / shortest-path / export live here ONCE; a
# concrete graph only supplies its nodes and edges.
#
# Per-edge `directed` distinguishes a directional relationship (a solve step, a
# `spectral_origin`) from a symmetric one (a Fourier pair, a co-constraint):
# traversal follows a directed edge `from → to` only, a symmetric edge both
# ways.  Stdlib-only; no rendering here (a view package consumes the data).

"""
    TypedEdge{N}(kind, from, to, detail, directed=true)

One typed edge of a [`KnowledgeGraph`](@ref): a `kind`-labeled connection from
node `from` to node `to` (both of type `N`), a human-readable `detail`, and
whether it is `directed` (traversed `from → to` only) or symmetric (both ways).
"""
struct TypedEdge{N}
    kind::Symbol
    from::N
    to::N
    detail::String
    directed::Bool
end
function TypedEdge(
    kind::Symbol, from::N, to::N, detail::AbstractString, directed::Bool=true
) where {N}
    return TypedEdge{N}(kind, from, to, String(detail), directed)
end

function Base.show(io::IO, e::TypedEdge)
    arrow = e.directed ? " → " : " — "
    return print(
        io, "TypedEdge(:", e.kind, ", ", e.from, arrow, e.to, ", \"", e.detail, "\")"
    )
end

"""
    KnowledgeGraph{N}(edges)

A knowledge graph over nodes of type `N`: a list of [`TypedEdge`](@ref)s.  The
generic parent the quantity / derivation / model graphs instantiate; iterate it
to get its edges (`for e in g`, `length(g)`, `collect(g)`).  Query it with
[`graph_nodes`](@ref), [`graph_neighbors`](@ref), [`graph_reachable`](@ref),
[`graph_shortest_path`](@ref); export it with [`graph_jsonl`](@ref).
"""
struct KnowledgeGraph{N}
    edges::Vector{TypedEdge{N}}
end

Base.iterate(g::KnowledgeGraph, s...) = iterate(g.edges, s...)
Base.length(g::KnowledgeGraph) = length(g.edges)
Base.eltype(::Type{KnowledgeGraph{N}}) where {N} = TypedEdge{N}
Base.isempty(g::KnowledgeGraph) = isempty(g.edges)

"""
    graph_edges(g) -> Vector{TypedEdge}

The edges of `g` (the backing vector).
"""
graph_edges(g::KnowledgeGraph) = g.edges

"""
    graph_nodes(g) -> Set

Every node that appears as an endpoint of some edge of `g`.
"""
function graph_nodes(g::KnowledgeGraph{N}) where {N}
    s = Set{N}()
    for e in g.edges
        push!(s, e.from)
        push!(s, e.to)
    end
    return s
end

"""
    graph_neighbors(g, node) -> Vector{TypedEdge}

Every edge of `g` incident on `node` (as `from` or `to`) — its neighborhood,
direction ignored.
"""
function graph_neighbors(g::KnowledgeGraph{N}, node) where {N}
    return filter(e -> e.from === node || e.to === node, g.edges)
end

# adjacency respecting per-edge direction: from a node, the (edge, other) steps
# it permits — a directed edge only `from → to`, a symmetric edge both ways.
function _adjacency(g::KnowledgeGraph{N}) where {N}
    adj = Dict{N,Vector{Tuple{TypedEdge{N},N}}}()
    for e in g.edges
        push!(get!(adj, e.from, Tuple{TypedEdge{N},N}[]), (e, e.to))
        e.directed || push!(get!(adj, e.to, Tuple{TypedEdge{N},N}[]), (e, e.from))
    end
    return adj
end

"""
    graph_reachable(g, start) -> Set

The set of nodes reachable from `start` following edge direction (directed
edges `from → to`, symmetric edges both ways); includes `start`.
"""
function graph_reachable(g::KnowledgeGraph{N}, start) where {N}
    adj = _adjacency(g)
    seen = Set{N}((start,))
    frontier = N[start]
    while !isempty(frontier)
        nxt = N[]
        for u in frontier, (_, v) in get(adj, u, Tuple{TypedEdge{N},N}[])
            v in seen && continue
            push!(seen, v)
            push!(nxt, v)
        end
        frontier = nxt
    end
    return seen
end

"""
    graph_shortest_path(g, a, b) -> Union{Vector{TypedEdge},Nothing}

A shortest path from node `a` to node `b` as the sequence of edges connecting
them (respecting per-edge direction); `nothing` if `b` is unreachable from `a`,
and the empty vector if `a === b`.
"""
function graph_shortest_path(g::KnowledgeGraph{N}, a, b) where {N}
    a === b && return TypedEdge{N}[]
    adj = _adjacency(g)
    prev = Dict{N,Tuple{TypedEdge{N},N}}()
    seen = Set{N}((a,))
    frontier = N[a]
    while !isempty(frontier)
        nxt = N[]
        for u in frontier, (e, v) in get(adj, u, Tuple{TypedEdge{N},N}[])
            v in seen && continue
            push!(seen, v)
            prev[v] = (e, u)
            v === b && return _rebuild_path(prev, a, b)
            push!(nxt, v)
        end
        frontier = nxt
    end
    return nothing
end

function _rebuild_path(prev::Dict{N,Tuple{TypedEdge{N},N}}, a, b) where {N}
    path = TypedEdge{N}[]
    node = b
    while node !== a
        e, u = prev[node]
        pushfirst!(path, e)
        node = u
    end
    return path
end

# hand-rolled JSON string (stdlib-only, matching QAtlas's *_jsonl idiom)
_json_str(s) = string('"', replace(string(s), '\\' => "\\\\", '"' => "\\\""), '"')

"""
    graph_jsonl([io=stdout], g; nodelabel=string) -> nothing

Stream `g` as JSONL for a network/graph view.  The first line is a summary
`{"nodes":N,"edges":M}`; each following line is one edge object
`{"kind":…,"from":…,"to":…,"detail":…,"directed":…}`, with node ids rendered
by `nodelabel` (default `string`; quantity graphs pass `nameof`).  No JSON
dependency — matching QAtlas's graph export so one consumer renders models ⊕
quantities ⊕ derivations.
"""
function graph_jsonl(io::IO, g::KnowledgeGraph; nodelabel=string)
    print(io, "{\"nodes\":", length(graph_nodes(g)), ",\"edges\":", length(g), "}\n")
    for e in g.edges
        print(
            io,
            "{\"kind\":",
            _json_str(e.kind),
            ",\"from\":",
            _json_str(nodelabel(e.from)),
            ",\"to\":",
            _json_str(nodelabel(e.to)),
            ",\"detail\":",
            _json_str(e.detail),
            ",\"directed\":",
            e.directed,
            "}\n",
        )
    end
    return nothing
end
graph_jsonl(g::KnowledgeGraph; kw...) = graph_jsonl(stdout, g; kw...)

export TypedEdge,
    KnowledgeGraph,
    graph_edges,
    graph_nodes,
    graph_neighbors,
    graph_reachable,
    graph_shortest_path,
    graph_jsonl
