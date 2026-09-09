using AbstractQAtlas
using AbstractQAtlas: check, conjugate_field, derivative_edge, index_spaces, tensor_rank
using Test

@testset "VectorPotential carries its dimension in the type" begin
    @test dimension(VectorPotential(0.3)) == 1
    @test dimension(VectorPotential(0.1, 0.2)) == 2
    @test dimension(VectorPotential(0.1, 0.2, 0.3)) == 3
    @test VectorPotential(0.3)[1] == 0.3
    @test length(VectorPotential(0.1, 0.2)) == 2
    @test collect(VectorPotential(0.1, 0.2)) == [0.1, 0.2]
    @test VectorPotential(1, 0.5) isa VectorPotential{2,Float64}
end

@testset "N is the embedding dimension, not the site set's" begin
    # A site set can have any Hausdorff dimension; the DISPLACEMENTS between its sites are
    # still vectors of the space it is drawn in, and that is what `A` contracts with. Three
    # displacements of a gasket-shaped set in the plane: the phases follow them and know
    # nothing about log3/log2.
    A = VectorPotential(0.3, -0.2)
    @test dimension(A) == 2
    for d in ((1.0, 0.0), (0.5, sqrt(3) / 2), (-0.5, sqrt(3) / 2))
        @test peierls_phase(A, d) ≈ 0.3 * d[1] - 0.2 * d[2]
    end
    # A cut-and-project chain lives on a line even though its hyperspace is 2D: its optical
    # `A` is one-dimensional, and a 2-vector cannot contract with its displacements.
    @test dimension(VectorPotential(0.3)) == 1
    @test_throws MethodError peierls_phase(VectorPotential(0.3, 0.0), (0.5,))
end

@testset "peierls_phase is A · d, and the dimensions must match" begin
    @test peierls_phase(VectorPotential(0.4), (0.5,)) ≈ 0.2
    @test peierls_phase(VectorPotential(0.1, 0.2), (1.0, 2.0)) ≈ 0.5
    @test_throws MethodError peierls_phase(VectorPotential(0.3), (0.5, 0.5))
    @test_throws MethodError peierls_phase(VectorPotential(0.1, 0.2), (0.5,))
    A = VectorPotential(0.3, -0.2)
    @test peierls_phase(A, (2.0, 4.0)) ≈ 2 * peierls_phase(A, (1.0, 2.0))
    @test peierls_phase(VectorPotential(0.6, -0.4), (1.0, 2.0)) ≈
        2 * peierls_phase(A, (1.0, 2.0))
    # The length unit lives in the displacement, which is the model's to supply: the same
    # field over half the distance is half the phase, and this layer names neither.
    @test peierls_phase(VectorPotential(0.4), (0.5,)) ≈
        peierls_phase(VectorPotential(0.4), (1.0,)) / 2
end

@testset "the current, the field and the potential are one genealogy edge" begin
    # Walked through the SAME accessors the thermodynamic quantities use, so a change to
    # either side breaks this — not "both happen to be 2".
    e = derivative_edge(ElectricCurrent)
    @test e.parent === Energy
    @test e.field === VectorPotentialField
    @test conjugate_field(ElectricCurrent()) === VectorPotentialField()
    @test VectorPotentialField() isa AbstractField
    # The magnetisation edge is the shape this copies, and it is NOT the same edge: `A`
    # couples to the hopping, so the parent is the Hamiltonian and not the free energy.
    @test derivative_edge(Magnetization(:z)).field === MagneticField
    @test derivative_edge(Magnetization(:z)).parent !== e.parent
end

@testset "the sign is the relation's, where every other signed derivative lives" begin
    # `j = -∂H/∂A`, stated as a relation exactly like `M = -∂F/∂h`, so it is checkable by
    # the same machinery rather than by a bare negation somewhere in this file.
    @test check(ElectricCurrentResponse(); j=-2.5, dH_dA=2.5)
    @test !check(ElectricCurrentResponse(); j=2.5, dH_dA=2.5)
    # A driven current oscillates about zero, so the sign is exactly what `|j|` cannot see:
    # both of the above have the same magnitude.
    @test abs(-2.5) == abs(2.5)
end

@testset "VectorPotential carries the same index traits as the quantities" begin
    for n in (1, 2, 3)
        A = VectorPotential(ntuple(i -> 0.1i, n))
        @test tensor_rank(typeof(A)) == n
        @test index_spaces(typeof(A)) == ntuple(_ -> SpatialDirection(), n)
        @test dimension(A) == tensor_rank(typeof(A))
    end
    @test tensor_rank(ElectricCurrent) == 1
    @test index_spaces(ElectricCurrent) == (SpatialDirection(),)
    @test index_spaces(typeof(VectorPotential(0.1, 0.2)))[1] ===
        index_spaces(ElectricCurrent)[1]
end

@testset "Region is the support, and is not a displacement" begin
    R = Region((1, 1), (1, 2))
    @test length(R) == 2
    @test_throws MethodError peierls_phase(VectorPotential(0.1, 0.2), R)
    @test peierls_phase(VectorPotential(0.1, 0.2), (1, 2)) ≈ 0.5
end

@testset "peierls_current refuses rather than answering a narrower question" begin
    A1 = VectorPotential(0.3)
    @test_throws ErrorException invoke(
        peierls_current, Tuple{Any,VectorPotential}, a -> a[1]^2, A1
    )
    msg1 = try
        invoke(peierls_current, Tuple{Any,VectorPotential}, a -> a[1]^2, A1)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("ForwardDiff", msg1)
    msg2 = try
        peierls_current(a -> 1.0, VectorPotential(0.1, 0.2))
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("2-dimensional", msg2)
    @test occursin("reverse-mode", msg2)
    @test !occursin("ForwardDiff", msg2)
end
