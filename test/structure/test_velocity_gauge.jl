using AbstractQAtlas
using Test

@testset "PeierlsConvention: the length unit" begin
    @test PeierlsConvention(1).bonds_per_cell == 1
    @test PeierlsConvention(2).bonds_per_cell == 2
    @test_throws ArgumentError PeierlsConvention(0)
    @test_throws ArgumentError PeierlsConvention(-1)
end

@testset "peierls_phase divides by the unit" begin
    @test peierls_phase(PeierlsConvention(1), 0.4) == 0.4
    @test peierls_phase(PeierlsConvention(2), 0.4) == 0.2
    @test peierls_phase(PeierlsConvention(3), 0.6) ≈ 0.2
    # Zero field carries no phase in any unit, and the map is linear in A — the two
    # properties a consumer relies on when it factors the phase out of a gate.
    for n in (1, 2, 4)
        c = PeierlsConvention(n)
        @test peierls_phase(c, 0.0) == 0.0
        @test peierls_phase(c, 2.0) ≈ 2 * peierls_phase(c, 1.0)
    end
    # The units are NOT interchangeable: this is the whole reason the type exists.
    @test peierls_phase(PeierlsConvention(1), 0.4) !=
        peierls_phase(PeierlsConvention(2), 0.4)
end

@testset "current_from_hamiltonian_derivative is one negation" begin
    @test current_from_hamiltonian_derivative(2.5) == -2.5
    @test current_from_hamiltonian_derivative(-2.5) == 2.5
    @test current_from_hamiltonian_derivative([1.0, -2.0]) == [-1.0, 2.0]
    @test current_from_hamiltonian_derivative([1.0 0.0; 0.0 2.0]) == [-1.0 0.0; 0.0 -2.0]
    # Not the identity, which is what the slip degrades to and what `|J|` cannot see.
    @test current_from_hamiltonian_derivative(2.5) != 2.5
    # An oscillating derivative comes back oscillating with the opposite phase, which is the
    # shape the sign error actually takes on a driven trace.
    d = [sin(t) for t in range(0, 2π; length=9)]
    @test current_from_hamiltonian_derivative(d) ≈ -d
    @test !isapprox(current_from_hamiltonian_derivative(d), d)
end

@testset "peierls_current refuses without a backend" begin
    # The seam errors rather than silently returning something: a missing AD backend must
    # not look like a zero current.
    if !isdefined(Main, :ForwardDiff)
        @test_throws ErrorException peierls_current(a -> a^2, 0.3)
    end
end
