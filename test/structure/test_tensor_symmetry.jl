# Intrinsic permutation symmetry of the nonlinear response tensors —
# the structural fact that makes the higher-order susceptibility a proper
# tensor: its field indices (with their frequencies) are interchangeable.

using AbstractQAtlas
using AbstractQAtlas:
    intrinsic_permutation_symmetric, canonical_component, permutation_equivalent

@testset "which quantities carry intrinsic permutation symmetry" begin
    @test intrinsic_permutation_symmetric(Susceptibility(:x, :y, :z))
    @test intrinsic_permutation_symmetric(DynamicalSusceptibility(:x, :y, :z))
    @test intrinsic_permutation_symmetric(Conductivity(:x, :y, :z))
    @test !intrinsic_permutation_symmetric(Energy())
    @test !intrinsic_permutation_symmetric(SpontaneousMagnetization())
    @test !intrinsic_permutation_symmetric(Magnetization(:x))   # rank-1, nothing to permute
end

@testset "canonical component sorts the FIELD indices (response index fixed)" begin
    # the response index (first) is fixed; field indices sort
    @test canonical_component(Susceptibility(:x, :z, :y)) === Susceptibility(:x, :y, :z)
    @test canonical_component(Susceptibility(:z, :z, :x)) === Susceptibility(:z, :x, :z)
    # linear (one field index): already canonical
    @test canonical_component(Susceptibility(:x, :y)) === Susceptibility(:x, :y)
    # dynamical and conductivity families too
    @test canonical_component(DynamicalSusceptibility(:x, :z, :y)) ===
        DynamicalSusceptibility(:x, :y, :z)
    @test canonical_component(Conductivity(:z, :y, :x)) === Conductivity(:z, :x, :y)
end

@testset "permutation equivalence" begin
    # χ⁽²⁾_{x;yz} = χ⁽²⁾_{x;zy}: same response index, permuted fields
    @test permutation_equivalent(Susceptibility(:x, :y, :z), Susceptibility(:x, :z, :y))
    # different RESPONSE index ⇒ not equivalent
    @test !permutation_equivalent(Susceptibility(:x, :y, :z), Susceptibility(:y, :x, :z))
    # different field content ⇒ not equivalent
    @test !permutation_equivalent(Susceptibility(:x, :y, :y), Susceptibility(:x, :y, :z))
    # reflexive on the linear component
    @test permutation_equivalent(Susceptibility(:x, :y), Susceptibility(:x, :y))
    # different families never equivalent
    @test !permutation_equivalent(Susceptibility(:x, :y, :z), Conductivity(:x, :y, :z))
    # non-symmetric quantities: not equivalent (no permutation freedom)
    @test !permutation_equivalent(Energy(), Energy())
end

@testset "field_permutation exposes the paired frequency permutation" begin
    using AbstractQAtlas: field_permutation, indices, frequency_arguments

    # the permutation that sorts the field indices IS the one the symmetry
    # applies to the frequency arguments (pairs (βᵢ,ωᵢ))
    @test field_permutation(DynamicalSusceptibility(:x, :y, :z)) == (1, 2)  # already canonical
    @test field_permutation(DynamicalSusceptibility(:x, :z, :y)) == (2, 1)  # swap the two ω's
    @test field_permutation(Susceptibility(:x, :z, :y)) == (2, 1)

    # INVARIANT: applying π to a component's field indices yields the
    # canonical component's field indices — so π is the frequency-argument
    # permutation needed to compare χ against canonical_component(χ).
    for χ in (
        DynamicalSusceptibility(:x, :z, :y),
        DynamicalSusceptibility(:z, :y, :x, :w),
        Susceptibility(:x, :z, :y),
        Conductivity(:z, :y, :x),
        DynamicalConductivity(:x, :z, :y),
    )
        π = field_permutation(χ)
        fields = collect(indices(χ)[2:end])
        @test Tuple(fields[collect(π)]) == indices(canonical_component(χ))[2:end]
        # one π per field slot = one per frequency argument (for the dynamical case)
        @test length(π) == length(fields)
    end
    # the dynamical response has as many frequencies as field slots to permute
    χd = DynamicalSusceptibility(:x, :z, :y)
    @test length(field_permutation(χd)) == frequency_arguments(χd)

    # no intrinsic symmetry ⇒ no field permutation
    @test_throws ErrorException field_permutation(Energy())
end
