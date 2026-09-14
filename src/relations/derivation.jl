# relations/derivation.jl — the relation registry as a DIRECTED derivation
# graph, and a lazy path-finding solver over it.
#
# Every `@relation` is an exact identity among its variables, and `solve`
# turns it into a computation: given all-but-one variable, produce the last.
# Read that way the whole registry is a DIRECTED graph — a node per relation
# variable (`:f`, `:Z`, `:β`, `:χ`, …), and, for every relation, a directed
# hyperedge `{other variables} →[relation] output` for each variable the
# relation can be solved for.  From a set of KNOWN quantities you can then ask
# what else is reachable, and — lazily — find one route to a target and run it.
#
# Two honesty guarantees, because a chained result is weaker than a directly
# implemented one:
#
#   * `derive` NEVER fabricates an edge: it calls the real `solve`, and a step
#     whose relation is non-affine in the wanted variable simply throws and is
#     skipped (the generic solver refuses non-affine variables by design).  So
#     the reachable set / route reflects what is ACTUALLY computable, not a
#     structural over-estimate.
#   * with `debug=true` the value comes wrapped in a [`DerivationTrace`] naming
#     the exact route (which relation produced each intermediate, from which
#     inputs) and flagged `indirect` — an indirectly derived number carries its
#     provenance so it is never mistaken for a directly-verified one.

"""
    DerivationStep(relation, output, inputs)

One directed edge of the derivation graph: `relation` computes the variable
`output::Symbol` from the variables `inputs::Tuple{Vararg{Symbol}}` (its other
variables), via [`solve`](@ref)`(relation, Val(output); inputs...)`.
"""
struct DerivationStep
    relation::AbstractRelation
    output::Symbol
    inputs::Tuple{Vararg{Symbol}}
end

function Base.show(io::IO, s::DerivationStep)
    ins = join(string.(s.inputs), ", ")
    return print(io, nameof(typeof(s.relation)), ": {", ins, "} → :", s.output)
end

"""
    DerivationTrace

The meta-information returned by [`derive`](@ref)`(...; debug=true)`: the
`target` symbol, its computed `value`, the ordered `steps` that produced it,
and `indirect` — `false` only when the target was among the supplied knowns
(a direct value), `true` when it was derived through one or more relations.
The trace exists for SAFETY: an indirectly derived value is auditable by its
route rather than trusted blindly.
"""
struct DerivationTrace
    target::Symbol
    value::Any
    steps::Vector{DerivationStep}
    indirect::Bool
end

function Base.show(io::IO, t::DerivationTrace)
    tag = t.indirect ? "indirect" : "direct"
    println(io, "DerivationTrace(:", t.target, " = ", t.value, "  [", tag, "])")
    if isempty(t.steps)
        print(io, "  (given directly)")
    else
        for (i, s) in enumerate(t.steps)
            print(io, "  ", i, ". ", s)
            i < length(t.steps) && println(io)
        end
    end
end

const _DERIV_STEPS = Ref{Union{Nothing,Vector{DerivationStep}}}(nothing)

"""
    derivation_steps() -> Vector{DerivationStep}

Every candidate directed edge of the derivation graph: for each registered
relation and each of its variables, the edge that would compute that variable
from the others.  These are STRUCTURAL candidates — a relation appears as an
edge for a variable it is not affine in too; [`derive`](@ref) is the honest
evaluator that discovers, by actually calling [`solve`](@ref), which edges
fire.  Built once and cached.
"""
function derivation_steps()
    cached = _DERIV_STEPS[]
    cached === nothing || return cached
    steps = DerivationStep[]
    for rel in all_relations()
        # EQUALITIES ONLY: `solve` on an inequality returns the SATURATION
        # (the tight bound), not an equational derivation — letting that into
        # the graph would derive a bound and present it as a computed value.
        rel isa AbstractInequality && continue
        vs = variables(rel)
        for out in vs
            ins = Tuple(v for v in vs if v !== out)
            push!(steps, DerivationStep(rel, out, ins))
        end
    end
    _DERIV_STEPS[] = steps
    return steps
end

# Try to fire a step against a value dict; return the solved value or `nothing`
# (the relation is non-affine in `output`, or the solve is otherwise refused).
function _try_step(step::DerivationStep, known::AbstractDict)
    all(v -> haskey(known, v), step.inputs) || return nothing
    try
        return solve(
            step.relation, Val(step.output); (v => known[v] for v in step.inputs)...
        )
    catch
        return nothing
    end
end

# Forward-chaining closure: keep firing any step whose inputs are all known
# until nothing new is produced.  Records the ordered steps actually used.
# Stops early once `stop` (if given) becomes known.
function _forward_chain(known::Dict{Symbol,Any}, stop::Union{Symbol,Nothing})
    used = DerivationStep[]
    progress = true
    while progress && !(stop !== nothing && haskey(known, stop))
        progress = false
        for step in derivation_steps()
            haskey(known, step.output) && continue
            v = _try_step(step, known)
            v === nothing && continue
            known[step.output] = v
            push!(used, step)
            progress = true
            stop !== nothing && haskey(known, stop) && break
        end
    end
    return used
end

# Keep only the steps that actually contribute to `target`, in dependency order.
function _prune(used::Vector{DerivationStep}, target::Symbol, given::Set{Symbol})
    by_output = Dict(s.output => s for s in used)
    needed = DerivationStep[]
    seen = Set{Symbol}()
    function visit(sym)
        (sym in given || sym in seen) && return nothing
        push!(seen, sym)
        step = get(by_output, sym, nothing)
        step === nothing && return nothing
        for inp in step.inputs
            visit(inp)
        end
        return push!(needed, step)   # post-order ⇒ inputs precede their consumer
    end
    visit(target)
    return needed
end

"""
    derivable(; knowns...) -> Set{Symbol}

The set of quantity variables COMPUTABLE from the supplied known values — the
honest reachability of the derivation graph.  Each element is a variable that
can be obtained, directly or through a chain of [`solve`](@ref)s, from the
knowns; the knowns themselves are included.

```julia
derivable(; Z = 2.0, β = 1.0)          # ⊇ Set([:Z, :β, :f])  — F = −β⁻¹ln Z reachable
```

Honest, not structural: a symbol appears only if some relation is actually
affine-solvable for it along the way (non-affine steps are skipped, exactly as
[`derive`](@ref) would skip them).
"""
function derivable(; knowns...)
    known = Dict{Symbol,Any}(pairs(knowns))
    _forward_chain(known, nothing)
    return Set(keys(known))
end

"""
    derive(target::Symbol; debug=false, knowns...) -> value | DerivationTrace

Lazily derive `target` from the supplied known values by finding ONE route
through the derivation graph and running it.  Returns the computed value; with
`debug=true` returns a [`DerivationTrace`](@ref) instead — the value plus the
exact route (which relation produced each intermediate, from which inputs) and
whether it was `indirect`.

```julia
derive(:f; Z = 2.0, β = 1.0)                 # -0.6931…   (F = −β⁻¹ ln Z)
derive(:f; Z = 2.0, β = 1.0, debug = true)   # DerivationTrace: 1. FreeEnergyFromZ: {Z, β} → :f
```

The route is discovered by forward chaining with the REAL [`solve`](@ref), so
a step whose relation is non-affine in its output is skipped, not faked.
Throws if `target` is not reachable from the knowns.
"""
function derive(target::Symbol; debug::Bool=false, knowns...)
    known = Dict{Symbol,Any}(pairs(knowns))
    given = Set(keys(known))
    if haskey(known, target)
        return if debug
            DerivationTrace(target, known[target], DerivationStep[], false)
        else
            known[target]
        end
    end
    used = _forward_chain(known, target)
    haskey(known, target) || error(
        "derive: :$target is not reachable from $(collect(given)) — no chain of " *
        "affine-solvable relations connects them (some steps may need a supplied " *
        "derivative or a value the knowns do not provide).",
    )
    debug || return known[target]
    return DerivationTrace(target, known[target], _prune(used, target, given), true)
end

"""
    derivation_graph() -> KnowledgeGraph{Symbol}

The derivation graph as a [`KnowledgeGraph`](@ref) instance, for a network VIEW
and structural inspection: one DIRECTED edge `input →[relation] output` per
(step, input), the simple-graph projection of the [`derivation_steps`](@ref)
hyperedges.

!!! warning "Structural, not computational"
    Each `@relation`'s output needs ALL of its inputs, but the projection fans
    each hyperedge out to one edge per input — so [`graph_reachable`](@ref) on
    this graph OVER-APPROXIMATES computability (a single known input already
    "reaches" the output, and non-affine outputs are edges too).  For the
    HONEST "can I actually compute it" use [`derivable`](@ref) / [`derive`](@ref),
    which require every input and call the real `solve`.  Use this graph for
    rendering and structural connectivity only.

Edge `kind` is the relation's name.
"""
function derivation_graph()
    edges = TypedEdge{Symbol}[]
    for step in derivation_steps()
        rname = Symbol(nameof(typeof(step.relation)))
        for inp in step.inputs
            push!(
                edges,
                TypedEdge(
                    rname, inp, step.output, string(nameof(typeof(step.relation))), true
                ),
            )
        end
    end
    return KnowledgeGraph(edges)
end

# ─── Type-keyed derivation (collision-proof) ────────────────────────────
#
# The symbol-keyed derivation above keys on private formula LETTERS, which COLLIDE
# across relations (`:S` = thermal entropy / thermopower / structure factor / vN
# entropy): `derive(:Π; S = entropy, T)` silently fires `KelvinRelation`, reusing an
# entropy as a Seebeck coefficient (the demonstrated bug the whole type-keyed
# redesign exists to kill).  Keying the graph on quantity TYPES ([`VariableKey`])
# makes that structurally impossible — `ThermalEntropy` is not `Thermopower`, so the
# step never fires.  Operates over the type-keyed relations (those with identity
# slots); the supplied (untyped) slots are provided as `extras`, exactly as for the
# bag verbs.  See docs/design/type-keyed-interface.md §6.

"""
    TypedStep(relation, output::VariableKey, inputs::Vector{VariableKey})

One directed edge of the TYPE-keyed derivation graph: `relation` computes the
identity variable `output` (a [`VariableKey`](@ref)) from its other identity slots
`inputs` (plus its supplied slots, provided as `extras`), via the type-keyed
[`solve`](@ref).
"""
struct TypedStep
    relation::AbstractRelation
    output::VariableKey
    inputs::Vector{VariableKey}
end

const _TYPED_DERIV_STEPS = Ref{Union{Nothing,Vector{TypedStep}}}(nothing)

"""
    typed_derivation_steps() -> Vector{TypedStep}

Every candidate directed edge of the type-keyed derivation graph: for each
type-keyed relation and each of its identity slots, the edge computing that slot's
TYPE from the other identity slots.  Inequalities and symbol-only relations are
skipped (an inequality gives a saturation bound, not a derivation; a symbol-only
relation has no typed slots).  Built once and cached.
"""
function typed_derivation_steps()
    cached = _TYPED_DERIV_STEPS[]
    cached === nothing || return cached
    steps = TypedStep[]
    for rel in all_relations()
        rel isa AbstractInequality && continue
        idslots = [(s, k) for (s, k) in variable_slots(rel) if k !== nothing]
        isempty(idslots) && continue
        for (osym, okey) in idslots
            ins = [VariableKey(k) for (s, k) in idslots if s !== osym]
            push!(steps, TypedStep(rel, VariableKey(okey), ins))
        end
    end
    _TYPED_DERIV_STEPS[] = steps
    return steps
end

# Is a type available in `known` — aliasing InverseTemperature ⇄ Temperature (so a
# derived T is recognized as "known" when β was supplied, and vice versa, and neither
# is ever re-derived as a duplicate literal key that would trip `_slot_value`'s
# "both present" guard mid-chain).
_known(@nospecialize(ty::Type), known::Bag) = _slot_value(ty, known) !== nothing

# Reject a bag that supplies BOTH β and T up front (mirrors the bag verbs' guard),
# so the aliasing-aware checks below never encounter both literal keys at once.
function _check_one_temperature(bag::Bag)
    haskey(bag, VariableKey(InverseTemperature)) &&
        haskey(bag, VariableKey(Temperature)) &&
        error("bag has both InverseTemperature and Temperature — pass exactly one")
    return nothing
end

# Fire a typed step against the known bag + `extras`; the solved value or `nothing`
# (an input type is missing, or the relation is non-affine in the output / a
# supplied slot is absent — the real `solve` decides, so nothing is faked).  The
# whole attempt is guarded, so an unexpectedly-erroring step is skipped, never fatal.
function _try_typed_step(step::TypedStep, known::Bag, extras)
    try
        all(k -> _known(k.type, known), step.inputs) || return nothing
        return solve(step.relation, step.output.type, known; extras...)
    catch
        return nothing
    end
end

# Forward-chaining closure over VariableKey nodes: fire any step whose inputs are all
# known until nothing new is produced (stopping early once `stop` is known).  All
# "known" tests are aliasing-aware, so `known` never accumulates both β and T.
function _typed_chain!(known::Bag, extras, stop::Union{VariableKey,Nothing})
    used = TypedStep[]
    progress = true
    while progress && !(stop !== nothing && _known(stop.type, known))
        progress = false
        for step in typed_derivation_steps()
            _known(step.output.type, known) && continue
            v = _try_typed_step(step, known, extras)
            v === nothing && continue
            known[step.output] = v
            push!(used, step)
            progress = true
            stop !== nothing && _known(stop.type, known) && break
        end
    end
    return used
end

"""
    derivable(bag::Bag; extras...) -> Set{VariableKey}

Type-keyed [`derivable`](@ref): the quantity/field TYPES computable from the bag
`bag` (plus `extras` for supplied slots), through chains of the type-keyed
[`solve`](@ref).  Matched by TYPE, so no formula-symbol collision — a
`ThermalEntropy` in the bag can never pose as a `Thermopower`.
"""
function derivable(bag::Bag; extras...)
    known = copy(bag)
    _check_one_temperature(known)
    _typed_chain!(known, values(extras), nothing)
    return Set(keys(known))
end

"""
    derive(Q::Type, bag::Bag; extras...) -> value

Type-keyed [`derive`](@ref): the value of quantity/field TYPE `Q` computed from
`bag` (+ `extras`), or an error if unreachable.  Collision-proof — a relation fires
only when the ACTUAL quantity types it needs are present:

```julia
derive(PeltierCoefficient, bag(Thermopower => s, Temperature => t))     # t·s (Kelvin)
derive(PeltierCoefficient, bag(ThermalEntropy => s, Temperature => t))  # ERROR: unreachable
#   — KelvinRelation needs a Thermopower, not the entropy; the silent
#   entropy-as-Seebeck derivation the symbol-keyed graph allowed is impossible here.
```
"""
function derive(@nospecialize(Q::Type), bag::Bag; extras...)
    known = copy(bag)
    _check_one_temperature(known)
    v = _slot_value(Q, known)                         # aliasing-aware "already known"
    v === nothing || return something(v)
    _typed_chain!(known, values(extras), VariableKey(Q))
    v = _slot_value(Q, known)
    v === nothing && error(
        "derive: $(nameof(Q)) is not reachable from the bag by type-keyed relations " *
        "(a needed quantity type or supplied slot is absent).",
    )
    return something(v)
end

"""
    typed_derivation_graph() -> KnowledgeGraph{VariableKey}

The type-keyed derivation graph: one directed edge `input →[relation] output` per
(typed step, input), nodes are [`VariableKey`](@ref)s — the collision-proof
counterpart of [`derivation_graph`](@ref).  Structural (over-approximates, like its
symbol sibling; use [`derive`](@ref)`(Q, bag)` for honest reachability).
"""
function typed_derivation_graph()
    edges = TypedEdge{VariableKey}[]
    for step in typed_derivation_steps()
        rname = Symbol(nameof(typeof(step.relation)))
        for inp in step.inputs
            push!(edges, TypedEdge(rname, inp, step.output, string(rname), true))
        end
    end
    return KnowledgeGraph(edges)
end

# ─── Quantity-first navigation (a veneer over the typed graph) ───────────

"""
    reachable_quantities(q) -> Vector{Type}

The physical-quantity TYPES structurally reachable from `q` through the type-keyed
derivation graph — the quantity-first navigation over [`typed_derivation_graph`](@ref).
Accepts a quantity instance or its (concrete) type; the result is family-erased
(`Susceptibility{I}` → `Susceptibility`), de-duplicated, name-sorted, and includes
`q`'s own family (trivially reachable).

```julia
reachable_quantities(PartitionFunction())   # ⊇ [FreeEnergy, PartitionFunction] — F = −β⁻¹ ln Z
```

The dual of [`relations_constraining`](@ref) (a quantity → the laws it obeys); this
is a quantity → the other quantities its laws connect it to.

!!! warning "Structural, not computational"
    Reachability over [`typed_derivation_graph`](@ref) OVER-approximates what is
    actually computable (a hyperedge fans out to one edge per input, so a single
    known input already "reaches" the output; non-affine outputs are edges too).
    For the honest "can I compute it from these values" use [`derivable`](@ref)`(bag)`
    / [`derive`](@ref)`(Q, bag)`, which require every input and call the real `solve`.
"""
reachable_quantities(q::AbstractQuantity) = reachable_quantities(typeof(q))
function reachable_quantities(@nospecialize(Q::Type))
    fam = _family(Q)
    self = fam <: AbstractQuantity ? Type[fam] : Type[]
    g = typed_derivation_graph()
    start = VariableKey(Q)
    start in graph_nodes(g) || return self          # not a node ⇒ only its own family
    reached = graph_reachable(g, start)
    qs = unique!(Type[_family(k.type) for k in reached if k.type <: AbstractQuantity])
    return sort!(qs; by=T -> string(nameof(T)))
end

export DerivationStep,
    DerivationTrace,
    derivation_steps,
    derivation_graph,
    derivable,
    derive,
    reachable_quantities,
    TypedStep,
    typed_derivation_steps,
    typed_derivation_graph

# ─── Cross-checking the routes, rather than taking one ─────────────────────
#
# `derive` finds one route and runs it. The registry usually offers several, and
# a quantity two relations both reach is a claim they have to agree on: each is
# an exact identity, so on one consistent set of knowns the answers coincide or
# one of the identities is wrong. That is a check no single relation can perform
# on itself, and it is what a hand-written cross-check does one target at a time.

"""
    ConsistencyRow

One held-out variable and every route back to it: the `target`, a `Symbol` on the
name-keyed report and a [`VariableKey`](@ref) on the type-keyed one, the value it
was held out at, the `steps` that reproduced it, their `values`, the `spread` over
those values and the held-out one, and whether they `agree`.
"""
struct ConsistencyRow
    target::Union{Symbol,VariableKey}
    held_out::Any
    steps::Vector{DerivationStep}
    values::Vector{Any}
    spread::Float64
    agree::Bool
end
export ConsistencyRow

function Base.show(io::IO, r::ConsistencyRow)
    return print(
        io,
        r.agree ? "agree" : "DISAGREE",
        " :",
        r.target,
        " over ",
        length(r.steps),
        " route(s), spread ",
        r.spread,
    )
end

"""
    consistency_report(data::NamedTuple; atol = 0, rtol = 1e-8, domain = nothing,
                       exclude = ()) -> Vector{ConsistencyRow}

Hold out each variable of `data` in turn and solve for it by every relation that
reaches it from the rest, then report whether those answers agree with each other
and with the value held out.

The registry is a set of exact identities, so a variable two of them both reach
must come back the same both ways and equal to what was removed. Where it does
not, one of the identities is wrong, and the row names every route so the odd one
out is visible. Nothing is hand-picked: the routes come from
[`derivation_steps`](@ref), the enumeration [`derive`](@ref) walks, and a step
non-affine in its target or otherwise refused drops out rather than being
counted, so a route appears only if it computes.

One step deep, from the remaining knowns. Chaining would compare a derived number
against another derived number, where a disagreement no longer names the relation
that caused it.

What is checked is the ALGEBRA, not the semantics. The registry is one namespace
and the report assumes a shared symbol is a shared quantity; where it is not, the
routes disagree and the row names them, but the fault is in the question rather
than in the identities. Two ways that happens, both measured:

  * A name means different things in different relations. `S` is a ring's block
    in one and an open chain's end block in another, so one number cannot satisfy
    both. `d` is the same trap: the classical image's dimension in
    [`Josephson`](@ref), the chain's own in [`HarrisCriterion`](@ref).
  * A relation does not apply at the point the data describes. Classical 2D Ising
    exponents agree over `:scaling` until `z` is added, at which point
    [`QuantumHyperscaling`](@ref) joins the routes to `α` and four rows disagree.
    It was excluded before only for want of an input, never for want of
    applicability, and nothing here knows the difference.

So a disagreement is a place to look, not a verdict. `domain` and `exclude` are
how a caller states the scope the data belongs to: `domain` keeps one family,
`exclude` drops named relations, and either is preferable to reading a row whose
routes describe different physics.

Agreement is `spread <= max(atol, rtol * scale)` with `scale` the largest
magnitude present, so `rtol` reads as a relative tolerance on the answer.

```julia
ising2d = (; α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4, d=2//1)
all(r -> r.agree, consistency_report(ising2d; domain=:scaling))
```
"""
function consistency_report(
    data::NamedTuple; atol=0, rtol=1e-8, domain::Union{Nothing,Symbol}=nothing, exclude=()
)
    out = ConsistencyRow[]
    steps = derivation_steps()
    domain === nothing ||
        (steps = filter(st -> AbstractQAtlas.domain(st.relation) === domain, steps))
    isempty(exclude) ||
        (steps = filter(st -> !(nameof(typeof(st.relation)) in exclude), steps))
    for target in keys(data)
        known = Dict{Symbol,Any}(k => v for (k, v) in pairs(data) if k !== target)
        got = Tuple{DerivationStep,Any}[]
        for st in steps
            st.output === target || continue
            v = _try_step(st, known)
            v === nothing && continue
            push!(got, (st, v))
        end
        isempty(got) && continue
        vals = [v for (_, v) in got]
        all(v -> v isa Number, vals) && data[target] isa Number || continue
        everything = vcat(float.(real.(vals)), float(real(data[target])))
        spread = maximum(everything) - minimum(everything)
        tol = max(atol, rtol * maximum(abs, everything))
        push!(
            out,
            ConsistencyRow(
                target, data[target], [st for (st, _) in got], vals, spread, spread <= tol
            ),
        )
    end
    return out
end
export consistency_report

"""
    consistent(data::NamedTuple; atol = 0, rtol = 1e-8, domain = nothing,
               exclude = ()) -> Bool

Whether every row of [`consistency_report`](@ref) agrees. `true` when no variable
is reachable at all, which is vacuous rather than a pass; read the report when
that matters.
"""
function consistent(data::NamedTuple; atol=0, rtol=1e-8, domain=nothing, exclude=())
    return all(
        r -> r.agree,
        consistency_report(data; atol=atol, rtol=rtol, domain=domain, exclude=exclude),
    )
end
export consistent

"""
    consistency_report(b::Bag; atol = 0, rtol = 1e-8, domain = nothing, exclude = (),
                       extras...) -> Vector{ConsistencyRow}

The type-keyed cross-check: the same leave-one-out sweep as the `NamedTuple`
method, over [`typed_derivation_steps`](@ref) instead of the symbol graph.

Sound where the other is not. A `VariableKey` is a quantity type and a support,
so two relations meet at a node only when they are talking about the same thing,
and the collision the name-keyed report cannot see does not arise: `β` is the
order-parameter exponent in `Rushbrooke` and the inverse temperature in
`DetailedBalance`, one symbol and two quantities, while `InverseTemperature` is a
type of its own. The typed known-set is aliasing-aware too, so a bag never holds
`InverseTemperature` and `Temperature` at once.

Narrower for the same reason. Only the 109 of 167 relations carrying at least one
typed slot appear at all, and a relation is present only through those slots, so
the classical scaling identities, whose exponents are bare symbols, are absent
here and are exactly what the name-keyed report checks best. Run both: this one
for what it can see, that one where a caller can vouch that a shared name is a
shared quantity.

Neither knows applicability; see the `NamedTuple` method.
"""
function consistency_report(
    b::Bag; atol=0, rtol=1e-8, domain::Union{Nothing,Symbol}=nothing, exclude=(), extras...
)
    out = ConsistencyRow[]
    steps = typed_derivation_steps()
    domain === nothing ||
        (steps = filter(st -> AbstractQAtlas.domain(st.relation) === domain, steps))
    isempty(exclude) ||
        (steps = filter(st -> !(nameof(typeof(st.relation)) in exclude), steps))
    for target in sort!(collect(keys(b)); by=k -> string(k.type))
        held = b[target]
        held isa Number || continue
        known = delete!(copy(b), target)
        got = Tuple{TypedStep,Any}[]
        for st in steps
            st.output.type === target.type || continue
            v = _try_typed_step(st, known, extras)
            v === nothing && continue
            push!(got, (st, v))
        end
        isempty(got) && continue
        vals = [v for (_, v) in got]
        all(v -> v isa Number, vals) || continue
        everything = vcat(float.(real.(vals)), float(real(held)))
        spread = maximum(everything) - minimum(everything)
        push!(
            out,
            ConsistencyRow(
                target,
                held,
                DerivationStep[
                    DerivationStep(
                        st.relation,
                        nameof(st.output.type),
                        Tuple(nameof(i.type) for i in st.inputs),
                    ) for (st, _) in got
                ],
                vals,
                spread,
                spread <= max(atol, rtol * maximum(abs, everything)),
            ),
        )
    end
    return out
end

"""
    consistent(b::Bag; kwargs...) -> Bool

Whether every row of the type-keyed [`consistency_report`](@ref) agrees.
"""
consistent(b::Bag; kwargs...) = all(r -> r.agree, consistency_report(b; kwargs...))
