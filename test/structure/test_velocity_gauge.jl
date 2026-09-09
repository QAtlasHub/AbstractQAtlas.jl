using AbstractQAtlas
using AbstractQAtlas: index_spaces
using Test

@testset "VectorPotential carries its dimension in the type" begin
    @test dimension(VectorPotential(0.3)) == 1
    @test dimension(VectorPotential(0.1, 0.2)) == 2
    @test dimension(VectorPotential(0.1, 0.2, 0.3)) == 3
    @test VectorPotential(0.3)[1] == 0.3
    @test length(VectorPotential(0.1, 0.2)) == 2
    @test collect(VectorPotential(0.1, 0.2)) == [0.1, 0.2]
    # Mixed reals promote rather than erroring, so a caller mixing an Int and a Float is not
    # forced to convert at the call site.
    @test VectorPotential(1, 0.5) isa VectorPotential{2,Float64}
end

@testset "N is the embedding dimension, not the set's" begin
    # A site set can have any Hausdorff dimension; the DISPLACEMENTS between its sites are
    # still ordinary vectors of the space it is drawn in, and that is what `A` contracts
    # with. Three sites of a Sierpiński-gasket-shaped set in the plane, with a 2-vector `A`:
    # the phases follow the displacements and know nothing about log3/log2.
    A = VectorPotential(0.3, -0.2)
    @test dimension(A) == 2
    for d in ((1.0, 0.0), (0.5, sqrt(3) / 2), (-0.5, sqrt(3) / 2))
        @test peierls_phase(A, d) ≈ 0.3 * d[1] - 0.2 * d[2]
    end
    # A cut-and-project chain lives on a line even though its hyperspace is 2D: its optical
    # `A` is one-dimensional, and a 2-vector cannot be contracted with its displacements.
    @test dimension(VectorPotential(0.3)) == 1
    @test_throws MethodError peierls_phase(VectorPotential(0.3, 0.0), (0.5,))
end

@testset "peierls_phase is A · d, and the dimensions must match" begin
    @test peierls_phase(VectorPotential(0.4), (0.5,)) ≈ 0.2
    @test peierls_phase(VectorPotential(0.1, 0.2), (1.0, 2.0)) ≈ 0.5
    # A bare number could not express this: a 1D potential and a 2D displacement is a type
    # error, not a silently truncated dot product.
    @test_throws MethodError peierls_phase(VectorPotential(0.3), (0.5, 0.5))
    @test_throws MethodError peierls_phase(VectorPotential(0.1, 0.2), (0.5,))
    # Linear in A and in d, which is what a consumer relies on to factor the phase out.
    A = VectorPotential(0.3, -0.2)
    @test peierls_phase(A, (2.0, 4.0)) ≈ 2 * peierls_phase(A, (1.0, 2.0))
    @test peierls_phase(VectorPotential(0.6, -0.4), (1.0, 2.0)) ≈
        2 * peierls_phase(A, (1.0, 2.0))
end

@testset "PeierlsConvention is the one-dimensional shortcut" begin
    @test bond_displacement(PeierlsConvention(1)) == 1
    @test bond_displacement(PeierlsConvention(2)) == 0.5
    @test_throws ArgumentError PeierlsConvention(0)
    @test_throws ArgumentError PeierlsConvention(-1)
    # It is the general phase with that displacement, not a second definition.
    for n in (1, 2, 4), a in (0.0, 0.4, -1.1)
        c, A = PeierlsConvention(n), VectorPotential(a)
        @test peierls_phase(c, A) == peierls_phase(A, (bond_displacement(c),))
        @test peierls_phase(c, A) ≈ a / n
    end
    # The units are NOT interchangeable: this is why the type exists.
    A = VectorPotential(0.4)
    @test peierls_phase(PeierlsConvention(1), A) != peierls_phase(PeierlsConvention(2), A)
end

@testset "current_from_hamiltonian_derivative is one negation" begin
    @test current_from_hamiltonian_derivative(2.5) == -2.5
    @test current_from_hamiltonian_derivative(-2.5) == 2.5
    @test current_from_hamiltonian_derivative([1.0, -2.0]) == [-1.0, 2.0]
    @test current_from_hamiltonian_derivative([1.0 0.0; 0.0 2.0]) == [-1.0 0.0; 0.0 -2.0]
    # Not the identity, which is what the slip degrades to and what `|J|` cannot see.
    @test current_from_hamiltonian_derivative(2.5) != 2.5
    d = [sin(t) for t in range(0, 2π; length=9)]
    @test current_from_hamiltonian_derivative(d) ≈ -d
    @test !isapprox(current_from_hamiltonian_derivative(d), d)
end

@testset "peierls_current refuses rather than answering a narrower question" begin
    # A missing AD backend must not read as "no current". Reached through `invoke` on the
    # fallback's own signature, so the assertion does not depend on whether ForwardDiff
    # happens to be loaded — a shard that ran the extension's tests first would otherwise
    # skip this silently, which is what the first version of this test did.
    A1 = VectorPotential(0.3)
    @test_throws ErrorException invoke(
        peierls_current, Tuple{Any,VectorPotential{1,Float64}}, a -> a[1]^2, A1
    )
    msg1 = try
        invoke(peierls_current, Tuple{Any,VectorPotential{1,Float64}}, a -> a[1]^2, A1)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("ForwardDiff", msg1)

    # A higher dimension is refused even WITH a backend, because reducing it to one
    # component would answer a different question and return a plausible number.
    msg2 = try
        peierls_current(a -> 1.0, VectorPotential(0.1, 0.2))
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("2-dimensional", msg2)
    @test occursin("reverse-mode", msg2)
    @test !occursin("ForwardDiff", msg2)   # not the backend message: a different reason
end

@testset "VectorPotential lines up with the index and support layers" begin
    # `N` counts SpatialDirection slots — the same space the transport tensors are indexed
    # by. A conductivity with n+1 spatial indices is driven by n fields, and a field for it
    # has one component per direction that space ranges over.
    σ = DynamicalConductivity(:x, :y)
    @test all(sp -> sp isa SpatialDirection, index_spaces(typeof(σ)))
    @test dimension(VectorPotential(0.1, 0.2)) == 2
    # A `Region` is the SUPPORT — which sites — and is deliberately label-based, so it
    # cannot supply a displacement. The two answer different questions, and the type system
    # keeps them apart: a Region's sites are not a tuple of reals.
    R = Region((1, 1), (1, 2))
    @test length(R) == 2
    @test_throws MethodError peierls_phase(VectorPotential(0.1, 0.2), R)
    # …and a displacement is not a Region either, so neither substitutes for the other.
    @test peierls_phase(VectorPotential(0.1, 0.2), (1, 2)) ≈ 0.5
end
