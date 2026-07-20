# fetch_cached — memoization for the retrieval seam (core-functions.md pillar 5).
# A mock model whose fetch method COUNTS its invocations lets us observe that a
# repeated lookup is served from the cache, that distinct keys recompute, and
# that errors are not cached.

using AbstractQAtlas
using AbstractQAtlas: fetch, fetch_cached, clear_fetch_cache!
using Test

struct _CacheModel <: AbstractQAtlas.AbstractQAtlasModel end
const _NCALLS = Ref(0)
function AbstractQAtlas.fetch(::_CacheModel, ::FreeEnergy, ::OBC; kwargs...)
    _NCALLS[] += 1
    return 42.0
end

@testset "fetch_cached memoizes the pure (model, quantity, bc) lookup" begin
    clear_fetch_cache!()
    _NCALLS[] = 0
    m = _CacheModel()
    v1 = fetch_cached(m, FreeEnergy(), OBC())
    v2 = fetch_cached(m, FreeEnergy(), OBC())
    @test v1 == v2 == 42.0
    @test _NCALLS[] == 1                          # computed exactly once — 2nd call is cached
    # a different kwarg is a different key ⇒ a genuine recomputation
    fetch_cached(m, FreeEnergy(), OBC(); extra=1)
    @test _NCALLS[] == 2
    # clearing the cache forces recomputation
    clear_fetch_cache!()
    fetch_cached(m, FreeEnergy(), OBC())
    @test _NCALLS[] == 3
end

@testset "fetch_cached forwards to fetch and does NOT cache errors" begin
    clear_fetch_cache!()
    # no fetch method for this (model, quantity, bc) ⇒ the erroring fallback fires;
    # the failure propagates and leaves the cache untouched
    @test_throws ErrorException fetch_cached(_CacheModel(), FreeEnergy(), PBC())
    # a subsequent working lookup is unaffected (nothing poisoned the cache)
    _NCALLS[] = 0
    @test fetch_cached(_CacheModel(), FreeEnergy(), OBC()) == 42.0
    @test _NCALLS[] == 1
    clear_fetch_cache!()
end
