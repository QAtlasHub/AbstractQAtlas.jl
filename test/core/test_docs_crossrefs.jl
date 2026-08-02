# Documenter cross-references, checked in the test suite instead of in the docs
# build — because two distinct dangling-@ref failures reached CI during the
# @bound work, and both were invisible to a `grep docs/`:
#
#   1. A stale `[`@inequality`](@ref)` inside the `AbstractInequality` DOCSTRING.
#      @autodocs renders docstrings onto index.md, so the reference is part of the
#      docs even though nothing under `docs/` mentions the macro.
#   2. A new physics submodule (`UniversalBounds`) absent from the @autodocs
#      `Modules =` array.  Its bindings resolve perfectly at the REPL, so
#      `isdefined` says everything is fine; its docstrings simply are not
#      rendered, and every @ref to them fails the build.
#
# The docs build catches both — 15 minutes later, in a separate workflow, on a
# push.  This catches them in the same run as everything else.

using AbstractQAtlas
using Test
const _M = AbstractQAtlas

_docs_root() = joinpath(pkgdir(AbstractQAtlas), "docs")

# The modules @autodocs will actually render, parsed from the block itself so
# this cannot drift from what the docs build does.
function _rendered_modules()
    index = read(joinpath(_docs_root(), "src", "index.md"), String)
    m = match(r"```@autodocs\s*\nModules\s*=\s*\[(.*?)\]"s, index)
    m === nothing && return nothing
    mods = Module[]
    for name in split(m.captures[1], ",")
        name = strip(name)
        isempty(name) && continue
        push!(mods, Core.eval(Main, Meta.parse(name)))
    end
    return mods
end

@testset "every physics submodule is rendered by @autodocs" begin
    listed = _rendered_modules()
    @test listed !== nothing
    @test _M in listed
    # not vacuous: the block really does enumerate the submodules
    @test length(listed) > 5
    for m in _M._PHYSICS_MODULES
        # a submodule missing here documents nothing, silently
        @test m in listed
    end
end

# A binding counts as documented if some rendered module carries a docstring for
# it — either the module it is defined in, or the one that owns it and re-exports.
function _documented(sym::Symbol, listed)
    for m in listed
        isdefined(m, sym) || continue
        haskey(Base.Docs.meta(m), Base.Docs.Binding(m, sym)) && return true
        owner = try
            parentmodule(getfield(m, sym))
        catch
            nothing
        end
        owner isa Module &&
            owner in listed &&
            haskey(Base.Docs.meta(owner), Base.Docs.Binding(owner, sym)) &&
            return true
    end
    return false
end

@testset "every [`x`](@ref) target resolves AND renders" begin
    listed = _rendered_modules()
    root = pkgdir(AbstractQAtlas)
    files = String[]
    for (dir, _, fs) in walkdir(joinpath(root, "src")), f in fs
        endswith(f, ".jl") && push!(files, joinpath(dir, f))
    end
    for f in readdir(joinpath(root, "docs", "src"); join=true)
        endswith(f, ".md") && push!(files, f)
    end
    @test length(files) > 20            # not vacuous

    pat = r"\[`([^`]+)`\]\(@ref\)"
    unresolved = String[]
    undocumented = String[]
    seen = Set{String}()
    for f in files, m in eachmatch(pat, read(f, String))
        # strip a call/parameter suffix: `fetch(model, q, bc)` -> `fetch`
        name = strip(split(split(m.captures[1], '(')[1], '{')[1])
        (isempty(name) || name in seen) && continue
        push!(seen, name)
        sym = Symbol(name)
        if !isdefined(_M, sym)
            isdefined(Base, sym) || push!(unresolved, "$(basename(f)): $name")
            continue
        end
        _documented(sym, listed) || push!(undocumented, "$(basename(f)): $name")
    end
    @test length(seen) > 50             # not vacuous: the docs really are cross-linked
    isempty(unresolved) || @info "@ref targets that name nothing" unresolved
    @test isempty(unresolved)
    isempty(undocumented) || @info "@ref targets with no RENDERED docstring" undocumented
    @test isempty(undocumented)
end
