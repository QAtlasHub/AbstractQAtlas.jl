#!/usr/bin/env julia
#
# test/ci/plan_shards.jl <N> [timings.tsv]
#
# Print (to stdout) a GitHub Actions matrix JSON that splits the test CASES across N
# shards, greedily balanced by recorded per-case timing (a `case<TAB>seconds` TSV)
# when one is supplied, else by count.  A CASE is either "aqua" (the module's Aqua
# suite) or a "dir/test_*.jl" file — the SAME list `test/runtests.jl` enumerates; each
# shard receives its subset via the `CI_CASES` env var (space-separated).  Matches the
# contract of the reusable `QAtlasHub/.github` `julia-ci.yml`:
#   [{"sid": "1", "cases": "aqua core/test_types.jl …"}, {"sid": "2", "cases": "…"}, …]
# Stdlib-only (no package deps): it runs before the package is built.

const DIRS = ["core", "structure", "relations", "report", "ext"]

# every runnable case: "aqua" + each dir/test_*.jl, sorted deterministically.
function all_cases()
    testdir = normpath(joinpath(@__DIR__, ".."))     # the test/ directory
    cases = String["aqua"]
    for d in DIRS
        dp = joinpath(testdir, d)
        isdir(dp) || continue
        for f in
            sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(dp)))
            push!(cases, string(d, "/", f))
        end
    end
    return cases
end

# per-case weight: seconds from `timings.tsv` (`case<TAB>seconds`) if present, else 1.0.
function weights(cases, timings)
    w = Dict(c => 1.0 for c in cases)
    (timings === nothing || !isfile(timings)) && return w
    for line in eachline(timings)
        parts = split(strip(line), '\t')
        length(parts) >= 2 || continue
        t = tryparse(Float64, parts[2])
        (t === nothing || !haskey(w, parts[1])) && continue
        w[parts[1]] = t
    end
    return w
end

# greedy longest-processing-time bin packing into N shards (balances total weight).
function pack(cases, w, n)
    bins = [String[] for _ in 1:n]
    load = zeros(Float64, n)
    for c in sort(cases; by=c -> -w[c])              # heaviest case first
        i = argmin(load)
        push!(bins[i], c)
        load[i] += w[c]
    end
    return filter(!isempty, bins)                    # never emit an empty shard
end

_json(s) = '"' * replace(String(s), '\\' => "\\\\", '"' => "\\\"") * '"'

function main()
    n = parse(Int, ARGS[1])
    timings = length(ARGS) >= 2 ? ARGS[2] : nothing
    cases = all_cases()
    bins = pack(cases, weights(cases, timings), max(1, min(n, length(cases))))
    objs = String[]
    for (i, b) in enumerate(bins)
        push!(
            objs,
            string(
                "{",
                _json("sid"),
                ":",
                _json(string(i)),
                ",",
                _json("cases"),
                ":",
                _json(join(b, " ")),
                "}",
            ),
        )
    end
    return print("[", join(objs, ","), "]")
end

main()
