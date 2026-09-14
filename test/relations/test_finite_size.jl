# test/relations/test_finite_size.jl — the SIZE axis: SizeSupport / at_size and the
# finite-size sweep that hands the supplied-derivative relations their derivative.

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
end

@testset "an activated sweep returns ψ exactly, and a conventional reading of it fails" begin
    # ln(1/Δ_typ) = A L^ψ with ψ = 1/2, which is the infinite-randomness form.
    ψ = 0.5
    Δ(L) = exp(-1.3 * L^ψ)
    sizes = (16, 32, 64, 128, 256)
    b = bag(
        (at_size(Typical(MassGap()), L) => Δ(L) for L in sizes)..., ActivatedExponent => ψ
    )
    rows = finite_size_scaling_report(b)

    # Four consecutive pairs from five sizes, and the secant is the derivative
    # exactly on this form, not to within a fit tolerance.
    @test length(rows) == 4
    @test all(r -> r.relation isa ActivatedFiniteSizeScaling, rows)
    @test all(r -> r.pass, rows)
    @test maximum(abs(r.residual) for r in rows) < 1e-12
    @test all(r -> r.quantity === Typical{MassGap}, rows)
    @test first(rows).sizes == (16, 32)

    # The same numbers read as a power of L instead of a power of ln L must fail,
    # or the report is not telling the two fixed points apart.
    conv = finite_size_scaling_report(
        bag(
            (at_size(Typical(MassGap()), L) => Δ(L) for L in sizes)...,
            DynamicalExponent => 2.0,
        ),
    )
    @test !any(r -> r.pass, conv)

    # And a wrong ψ has to fail too, or passing says nothing.
    @test !any(
        r -> r.pass,
        finite_size_scaling_report(
            bag(
                (at_size(Typical(MassGap()), L) => Δ(L) for L in sizes)...,
                ActivatedExponent => 0.8,
            ),
        ),
    )
end

@testset "a conventional sweep is read the other way, and both arms can be asked" begin
    z = 2.0
    Ω(L) = 0.7 * L^(-z)
    sizes = (8, 16, 32, 64)
    b = bag((at_size(MassGap(), L) => Ω(L) for L in sizes)..., DynamicalExponent => z)
    rows = finite_size_scaling_report(b)
    @test length(rows) == 3
    @test all(r -> r.relation isa ConventionalFiniteSizeEnergy, rows)
    @test all(r -> r.pass, rows)

    # A bag carrying both exponents is checked against both laws on the same
    # numbers, which is the comparison a size sweep is usually run to settle.
    both = bag(
        (at_size(Typical(MassGap()), L) => Ω(L) for L in sizes)...,
        DynamicalExponent => z,
        ActivatedExponent => 0.5,
    )
    kinds = Dict{Symbol,Vector{Bool}}()
    for r in finite_size_scaling_report(both)
        push!(get!(kinds, nameof(typeof(r.relation)), Bool[]), r.pass)
    end
    @test all(kinds[:ConventionalFiniteSizeEnergy])
    @test !any(kinds[:ActivatedFiniteSizeScaling])
end

@testset "ψ may not be read off an average, or off a sweep that does not say" begin
    Δ(L) = exp(-1.3 * sqrt(L))
    avg = bag(
        (at_size(DisorderAveraged(MassGap()), L) => Δ(L) for L in (16, 32))...,
        ActivatedExponent => 0.5,
    )
    # The average of a gap at an infinite-randomness fixed point is set by the rare
    # regions, so a secant off it is not ψ. Refused rather than skipped: skipping
    # would leave a caller wondering where their rows went.
    @test_throws "set by the rare regions" finite_size_scaling_report(avg)

    bare = bag(
        (at_size(MassGap(), L) => Δ(L) for L in (16, 32))..., ActivatedExponent => 0.5
    )
    @test_throws "say which this is" finite_size_scaling_report(bare)

    # Without ψ the same averaged sweep is read conventionally, no refusal: the
    # guard is about the activated law, not about averages.
    @test !isempty(
        finite_size_scaling_report(
            bag(
                (at_size(DisorderAveraged(MassGap()), L) => Δ(L) for L in (16, 32))...,
                DynamicalExponent => 2.0,
            ),
        ),
    )
end

@testset "the sweep reports nothing rather than something wrong" begin
    Δ(L) = exp(-1.3 * sqrt(L))
    # One size is not a sweep, and no exponent is nothing to check against.
    @test isempty(
        finite_size_scaling_report(
            bag(at_size(Typical(MassGap()), 16) => Δ(16), ActivatedExponent => 0.5)
        ),
    )
    @test isempty(
        finite_size_scaling_report(
            bag((at_size(Typical(MassGap()), L) => Δ(L) for L in (16, 32))...)
        ),
    )

    # `ln[-ln O]` needs 0 < O < 1. Values outside it get no activated row rather
    # than a row built on a complex or infinite intermediate.
    big = bag(
        (at_size(Typical(MassGap()), L) => 1.0 + L for L in (16, 32))...,
        ActivatedExponent => 0.5,
    )
    @test isempty(finite_size_scaling_report(big))
    @test length(
        finite_size_scaling_report(
            bag(
                (at_size(Typical(MassGap()), L) => 1.0 + L for L in (16, 32))...,
                DynamicalExponent => 2.0,
            ),
        ),
    ) == 1
end
