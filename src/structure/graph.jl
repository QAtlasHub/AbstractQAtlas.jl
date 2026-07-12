# structure/graph.jl — the quantity-relationship graph, an instance of the
# generic KnowledgeGraph kernel (graph/kernel.jl).
#
# AbstractQAtlas grows several *separate* typed edges between physical
# quantities, each historically its own ad-hoc accessor:
#
#   * `derivative_edge`          — the response genealogy (χ = ∂M/∂h, …),
#   * `spectral_origin`          — the dynamical graph (A = −ImG/π, …),
#   * `fourier_conjugate_quantity` — conjugate pairs (S(q) ↔ ⟨SS⟩(r)),
#   * the relation ↔ quantity links (`quantities` / `relations_constraining`
#     from quantity_links.jl) — every universal law a quantity must obey.
#
# They share one shape — a typed edge between quantity kinds — so they fold into
# one `KnowledgeGraph{Type}`: `related_quantities(q)` gives a quantity's
# neighborhood, and the generic kernel answers "how are A and B related?"
# (`quantity_path` ⟶ `graph_shortest_path`) and exports the network
# (`quantity_graph_jsonl` ⟶ `graph_jsonl`).  Nodes are quantity FAMILIES (the
# index-erased UnionAll — `Susceptibility`, not `Susceptibility{(:z,:z)}`) so
# the graph is finite; the concrete index still RESOLVES an edge (χ⁽²⁾ ⟶ χ⁽¹⁾ ⟶
# M) before the endpoints are normalized.  The structural edges are stored
# `directed=false` (relatedness is symmetric); their `from → to` orientation is
# kept for the view.

"""
    QuantityEdge

A [`TypedEdge`](@ref) of the quantity-relationship graph (an alias for
`TypedEdge{Type}`, nodes are quantity families).  Its `kind` is one of

- `:derivative` — `to` is the potential/quantity `from` is a field-derivative
  of (the response genealogy, [`derivative_edge`](@ref)); `detail` names the
  field (`"∂/∂MagneticField"`, …).
- `:spectral` — `from` is obtained from `to` by the dynamical-graph operation
  in `detail` ([`spectral_origin`](@ref)'s `via`: `dyson`, `neg_im_over_pi`, …).
- `:fourier` — `from` and `to` are Fourier conjugates
  ([`fourier_conjugate_quantity`](@ref)); `detail == "fourier"`.
- `:law` — `from` and `to` are co-constrained by a universal relation (they
  appear together in some [`quantities`](@ref)`(rel)`); `detail` names the
  relation (`"SusceptibilityFDT"`, …).
"""
const QuantityEdge = TypedEdge{Type}

# The index-erased FAMILY of a quantity type: `Susceptibility{(:z,:z)}` and the
# bare `Susceptibility` both normalize to the `Susceptibility` UnionAll, so the
# structural graph has one node per quantity kind.
_family(::Type{T}) where {T} = Base.typename(T).wrapper
_family(q::AbstractQuantity) = _family(typeof(q))

# A short human label for a field type in a derivative edge.
_field_label(::Type{F}) where {F} = string(nameof(F))

# Instantiate a representative of a (possibly parametric) quantity family, or
# `nothing` if none of the standard arities construct — used to probe the
# per-type edge accessors and to seed graph traversal.
function _rep_quantity(T::Type)
    for args in ((), (:x,), (:x, :y))
        try
            return T(args...)
        catch
        end
    end
    return nothing
end

"""
    related_quantities(q) -> Vector{QuantityEdge}

The graph neighborhood of quantity `q` (a quantity instance or type): every
[`QuantityEdge`](@ref) it participates in across ALL edge kinds —
its response-genealogy parent (`:derivative`), its dynamical origin
(`:spectral`), its Fourier conjugate (`:fourier`), and every quantity it is
co-constrained with by a universal law (`:law`).  Endpoints are quantity
families (the index-erased UnionAll).

```julia
related_quantities(Susceptibility(:z, :z))
# TypedEdge(:derivative, Susceptibility — Magnetization, "∂/∂MagneticField")
# TypedEdge(:law, Susceptibility — Magnetization, "SusceptibilityFDT")
# …
```
"""
related_quantities(q::AbstractQuantity) = _related_quantities(q, typeof(q))
function related_quantities(T::Type)
    # A concrete type dispatches the per-type accessors directly; a bare family
    # (UnionAll) would fall through to their `<:AbstractQuantity` defaults and
    # miss its index-parametric edge, so instantiate a representative first.
    isconcretetype(T) && return _related_quantities(T, T)
    r = _rep_quantity(T)
    return r === nothing ? _related_quantities(T, T) : _related_quantities(r, typeof(r))
end

# `probe` is a quantity instance OR type accepted by the per-type accessors;
# `ct` is the concrete type used to resolve index-dependent edges.  Structural
# edges are symmetric (`directed=false`); orientation is kept in from/to.
function _related_quantities(probe, ct::Type)
    edges = QuantityEdge[]
    self = _family(ct)

    e = derivative_edge(probe)
    if e !== nothing
        push!(
            edges,
            QuantityEdge(
                :derivative, self, _family(e.parent), "∂/∂$(_field_label(e.field))", false
            ),
        )
    end

    o = spectral_origin(probe)
    if o !== nothing
        push!(edges, QuantityEdge(:spectral, self, _family(o.from), string(o.via), false))
    end

    c = fourier_conjugate_quantity(ct)
    if c !== nothing
        push!(edges, QuantityEdge(:fourier, self, _family(c), "fourier", false))
    end

    # co-constraint edges: for every relation naming this quantity, link it to
    # each OTHER quantity that relation names (a shared universal law).
    for rel in relations_constraining(probe)
        rname = string(nameof(typeof(rel)))
        for U in quantities(rel)
            fam = _family(U)
            fam === self && continue
            push!(edges, QuantityEdge(:law, self, fam, rname, false))
        end
    end
    return edges
end

# recursive concrete leaves of the AbstractQuantity hierarchy
function _quantity_leaves(T::Type=AbstractQuantity, acc::Vector{Type}=Type[])
    for S in subtypes(T)
        isabstracttype(S) ? _quantity_leaves(S, acc) : push!(acc, S)
    end
    return acc
end

const _GRAPH_CACHE = Ref{Union{Nothing,KnowledgeGraph{Type}}}(nothing)

"""
    quantity_graph() -> KnowledgeGraph{Type}

The entire quantity-relationship graph as a [`KnowledgeGraph`](@ref) — every
`:derivative`, `:spectral`, `:fourier` and `:law` edge among the constructible
quantity families, deduplicated.  Built once from the type hierarchy and
cached; query it with the generic kernel ([`graph_neighbors`](@ref),
[`graph_shortest_path`](@ref), …) or the quantity-specific wrappers below.
"""
function quantity_graph()
    cached = _GRAPH_CACHE[]
    cached === nothing || return cached
    seen = Set{Tuple{Symbol,Type,Type,String}}()
    edges = QuantityEdge[]
    for T in _quantity_leaves()
        q = _rep_quantity(T)
        q === nothing && continue
        for e in related_quantities(q)
            # canonical unordered key so an edge and its mirror collapse
            key = (e.kind, e.from, e.to, e.detail)
            rkey = (e.kind, e.to, e.from, e.detail)
            (key in seen || rkey in seen) && continue
            push!(seen, key)
            push!(edges, e)
        end
    end
    g = KnowledgeGraph(edges)
    _GRAPH_CACHE[] = g
    return g
end

"""
    quantity_neighbors(fam) -> Vector{QuantityEdge}

Every edge of [`quantity_graph`](@ref) incident on quantity family `fam` (as
`from` OR `to`) — its neighborhood.  Unlike [`related_quantities`](@ref) this
also surfaces edges naming `fam` as their `to` endpoint (e.g. the quantities
whose derivative is `fam`).  A thin wrapper over [`graph_neighbors`](@ref).
"""
quantity_neighbors(fam::Type) = graph_neighbors(quantity_graph(), _family(fam))

"""
    quantity_path(a, b) -> Union{Vector{QuantityEdge},Nothing}

A shortest path in the quantity-relationship graph from quantity `a` to `b`
(instances, families, or types) — the machine answer to "how are `a` and `b`
related?".  `nothing` if they are in different components, the empty vector if
`a` and `b` are the same family.  A thin wrapper over
[`graph_shortest_path`](@ref) (structural edges are symmetric, so the search is
undirected).

```julia
quantity_path(SpecificHeat(), Magnetization(:z))
# SpecificHeat — Energy — FreeEnergy — Magnetization   (through the common root)
```
"""
function quantity_path(a, b)
    src = _family(a isa Type ? a : typeof(a))
    dst = _family(b isa Type ? b : typeof(b))
    return graph_shortest_path(quantity_graph(), src, dst)
end

"""
    quantity_graph_jsonl([io=stdout]) -> nothing

Stream the whole [`quantity_graph`](@ref) as JSONL for a network/graph view
(node ids are the family type names).  A thin wrapper over [`graph_jsonl`](@ref).
"""
quantity_graph_jsonl(io::IO=stdout) = graph_jsonl(io, quantity_graph(); nodelabel=nameof)

export QuantityEdge,
    related_quantities,
    quantity_graph,
    quantity_neighbors,
    quantity_path,
    quantity_graph_jsonl
