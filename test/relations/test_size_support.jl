# test/relations/test_size_support.jl — the SIZE axis of the support layer.

using AbstractQAtlas
using Test

@testset "SizeSupport is what makes a finite-size sweep expressible" begin
    # Without a support the sizes collide, and the bag says so rather than keeping
    # the last one. This is the whole reason the size axis exists.
    @test_throws "duplicate key" bag(Typical(MassGap()) => 1e-2, Typical(MassGap()) => 1e-4)
    b = bag(
        at_size(Typical(MassGap()), 16) => 1e-2, at_size(Typical(MassGap()), 32) => 1e-4
    )
    @test length(b) == 2
    @test b[at_size(Typical(MassGap()), 32)] == 1e-4

    # Value-keyed, like every other support: the same size is the same slot.
    @test at_size(Typical(MassGap()), 16) == at_size(Typical(MassGap()), 16)
    @test at_size(Typical(MassGap()), 16) != at_size(Typical(MassGap()), 32)
    @test at_size(MassGap, 16) == at_size(MassGap(), 16)

    # A key carries one support, so a quantity that already needs its own cannot also
    # be keyed by size: `typeof` erases the order, and two Rényi entropies at
    # different sizes would otherwise fuse into a sweep that never existed.
    @test_throws "cannot" at_size(RenyiEntropy(2.0), 16)
end
