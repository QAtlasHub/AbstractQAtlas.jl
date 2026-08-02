# Citation discipline — docs/design/core-functions.md, pillar 3.
#
# Every relations source file that DECLARES a physics relation must carry a
# citation: a `[key](@cite)` to the doiget-verified bibliography
# (docs/references.bib, checked by the `citations` CI gate) OR an honest inline
# `(YYYY)` textbook / historical reference, in the file header's "References:"
# block or a docstring.  A file that ships a law with NO provenance fails here —
# a verify-engine that cannot cite its own laws cannot be trusted (ties the
# never-fabricate-citations rule).  Per-relation `[key](@cite)` is the
# aspirational form; the enforced minimum is per-file, matching the established
# header-"References:" convention (e.g. scaling.jl).
#
# `interface.jl` is exempt: it DEFINES the @relation/@bound macros, and its
# `@relation :scaling Rushbrooke(...)` lines are docstring EXAMPLES (the real
# Rushbrooke is declared in scaling.jl), not registered relations.

using AbstractQAtlas
using Test

const _RELATIONS_DIR = joinpath(pkgdir(AbstractQAtlas), "src", "relations")
const _CITE_EXEMPT = ("interface.jl",)   # macro-definition home (its @relation lines are doc examples)

_declares_relation(src) = occursin(r"@(relation|bound)\s+:", src)
_has_citation(src) = occursin("(@cite)", src) || occursin(r"\((1[89]\d\d|20\d\d)\)", src)

@testset "every relations file cites its sources (core-functions.md pillar 3)" begin
    files = sort(filter(f -> endswith(f, ".jl"), readdir(_RELATIONS_DIR; join=true)))
    @test !isempty(files)
    for f in files
        basename(f) in _CITE_EXEMPT && continue
        src = read(f, String)
        _declares_relation(src) || continue         # only files that declare a relation
        cited = _has_citation(src)
        cited || @warn "relations file declares a law with NO citation" file = basename(f)
        @test cited
    end
end
