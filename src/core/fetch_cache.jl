# core/fetch_cache.jl — memoization for the `fetch` retrieval seam.
#
# A reference value is a PURE function of its `(model, quantity, bc)` key (the
# `fetch` contract), and the oracle behind it (ED / DMRG / TPQ) is expensive.
# `fetch_cached` wraps `fetch` with a stdlib `Dict` cache so a repeated lookup is
# free — the retrieval-efficiency seam (core-functions.md, pillar 5).  The cache
# lives HERE (the parent owns the seam); the expensive `fetch` METHODS live
# downstream (QAtlas).  Mirrors the build-once derivation-graph cache, extended
# to the value seam.  Thread-safe: the value is computed OUTSIDE the lock, so
# distinct keys compute concurrently and only the store is serialized.

struct _FetchMiss end
const _FETCH_MISS = _FetchMiss()                # sentinel distinct from any real value
const _FETCH_CACHE = Dict{Any,Any}()
const _FETCH_CACHE_LOCK = ReentrantLock()

"""
    fetch_cached(model, quantity, bc; kwargs...)

Memoizing wrapper over [`fetch`](@ref): the first call for a given
`(model, quantity, bc)` (plus `kwargs`) computes the value via `fetch` and stores
it; later identical calls return the stored value without recomputing.

Safe because a reference value is a pure function of its key (the `fetch`
contract), and the oracle behind it is expensive.  The key compares the arguments
**by value**, which is exactly right for the immutable model / quantity / BC
structs (a mutable model caches per object identity).  A `fetch` method with side
effects must not be cached.  Errors are **not** cached — a failed `fetch`
propagates and leaves the cache untouched.  Clear with [`clear_fetch_cache!`](@ref).

Thread-safe: the value is computed outside the lock so distinct keys compute
concurrently; a race on the same key keeps the first stored value.
"""
function fetch_cached(
    model::AbstractQAtlasModel, quantity::AbstractQuantity, bc::BoundaryCondition; kwargs...
)
    key = (model, quantity, bc, (; kwargs...))
    hit = Base.@lock _FETCH_CACHE_LOCK get(_FETCH_CACHE, key, _FETCH_MISS)
    hit === _FETCH_MISS || return hit
    value = fetch(model, quantity, bc; kwargs...)    # outside the lock (may be slow, may error)
    return Base.@lock _FETCH_CACHE_LOCK get!(_FETCH_CACHE, key, value)
end
export fetch_cached

"""
    clear_fetch_cache!()

Empty the [`fetch_cached`](@ref) memoization cache — e.g. after a model's
parameters change, invalidating stored values.  Returns `nothing`.
"""
function clear_fetch_cache!()
    Base.@lock _FETCH_CACHE_LOCK empty!(_FETCH_CACHE)
    return nothing
end
export clear_fetch_cache!
