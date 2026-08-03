ENV["GKSwstype"] = "100"

using AbstractQAtlas
using Test, Aqua
using TestShards

# Every `test_*.jl` under `test/`, in a deterministic order, each one its own shardable unit,
# plus Aqua. `@shard` shadows `include` inside the block, so a unit is whatever this loop
# includes — a new file, or a whole new directory, is picked up BY BEING ON DISK, rather than by
# being added to the `dirs` list and the `all_cases()` function that used to sit here and could
# both disagree with the tree.
#
# Two rules when adding to this, and they are the only two:
#
#   1. SHARED FIXTURES GO ABOVE THIS BLOCK. A helper included inside becomes a unit of its own,
#      lands on ONE shard, and every test file on the other shards that needed it fails.
#   2. ANYTHING THAT IS NOT A `test_*.jl` FILE MUST BE NAMED, as Aqua is below. The glob does not
#      error on what it does not match; it silently stops running it.
#
# A bare `Pkg.test()` with nothing set in the environment runs all of it, in this order. Run one
# shard locally with `TESTSHARDS_ID=s2 TESTSHARDS_N=4 julia --project -e 'using Pkg; Pkg.test()'`.
TestShards.@shard begin
    for (dir, _, files) in sort!(collect(walkdir(@__DIR__)))
        for f in sort(files)
            startswith(f, "test_") && endswith(f, ".jl") || continue
            include(joinpath(dir, f))
        end
    end
    # The module's own QA. Not a file, so `@unit` gives it a key of its own — it was already a
    # first-class case here (`all_cases()` began with "aqua"), and it stays one.
    TestShards.@unit "aqua" begin
        Aqua.test_all(AbstractQAtlas)
    end
end
