using AbstractQAtlas
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
