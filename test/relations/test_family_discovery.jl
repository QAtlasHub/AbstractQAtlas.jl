# test/relations/test_family_discovery.jl — §8a family-generic slot auto-discovery.
#
# A slot keyed on a bare parametric FAMILY (`Susceptibility`, a UnionAll) matches
# every concrete component of that family present in a bag; `relation_report` /
# `check_all` auto-instantiate one row per component. SusceptibilityPositivity is
# the prototype. See docs/design/type-keyed-interface.md §8a.

using AbstractQAtlas
using Test
const AQ = AbstractQAtlas

@testset "§8a: family-generic auto-discovery (SusceptibilityPositivity)" begin
    # keyed on the bare Susceptibility family; quantities still family-erase to it
    @test quantities(SusceptibilityPositivity()) == (Susceptibility,)
    @test AQ._is_family(only(variable_types(SusceptibilityPositivity())))

    # the symbol path is unchanged
    @test check(SusceptibilityPositivity(); χT=0.7)
    @test !check(SusceptibilityPositivity(); χT=-0.1, atol=1e-9)

    # auto-discovery: one bag, many components ⇒ one report row each
    b = bag(
        Susceptibility{(:x, :x)} => 1.0,
        Susceptibility{(:z, :z)} => -0.5,
        Susceptibility{(:y, :y)} => 2.0,
    )
    rep = relation_report(b; atol=1e-9)
    @test length(rep) == 3
    @test count(r -> r.pass, rep) == 2               # χ_xx, χ_yy ≥ 0
    @test count(r -> !r.pass, rep) == 1              # χ_zz < 0 — the violation is caught
    # every row names the concrete component it was instantiated on
    @test all(r -> r.subject isa VariableKey && r.subject.type <: Susceptibility, rep)
    @test Set(r.subject.type for r in rep) == Set((
        Susceptibility{(:x, :x)}, Susceptibility{(:z, :z)}, Susceptibility{(:y, :y)}
    ))
    # the whole bag fails, since a negative component exists (whichever it is)
    @test !check_all(b; atol=1e-9)

    # an explicit component is selected with `subject`
    @test check(SusceptibilityPositivity(), b; subject=Susceptibility{(:x, :x)})
    @test !check(SusceptibilityPositivity(), b; subject=Susceptibility{(:z, :z)}, atol=1e-9)
    # …and the bare bag form (ambiguous across components) is a loud error
    @test_throws ErrorException residual(
        SusceptibilityPositivity(), bag(Susceptibility{(:x, :x)} => 1.0)
    )
end

@testset "§8a: ordinary relations keep subject === nothing" begin
    gr = 1 / (0.2 + im * 0.1)
    rep = relation_report(
        bag(
            RetardedGreensFunction => gr,
            AdvancedGreensFunction => conj(gr),
            SpectralFunction => -imag(gr) / π,
        );
        atol=1e-12,
    )
    @test !isempty(rep)
    @test all(r -> r.subject === nothing, rep)       # concrete relations: no auto-instantiation
end

@testset "§8a review #86 regression: concrete-only + loud guards" begin
    # HIGH: solve rejects a bare FAMILY target (concrete-only; a family is ambiguous)
    @test_throws ErrorException solve(
        SusceptibilityPositivity(), Susceptibility, bag(SpecificHeat => 1.0)
    )
    # a family target is not derivable through the typed graph either
    @test_throws ErrorException derive(Susceptibility, bag(InverseTemperature => 2.0))

    # MEDIUM: a bag entry keyed on the bare family ITSELF is not a component…
    @test isempty(AQ._bag_components(Susceptibility, bag(Susceptibility => 3.0)))
    # …so it produces no bogus report row (subject must be a concrete component)
    @test isempty(relation_report(bag(Susceptibility => 3.0); atol=1e-9))

    # MEDIUM: `subject` on a NON-family-generic relation is a loud error, not ignored
    @test_throws ErrorException residual(
        SpecificHeatPositivity(), bag(SpecificHeat => 1.5); subject=Susceptibility{(:x, :x)}
    )
    @test_throws ErrorException check(
        SpecificHeatPositivity(), bag(SpecificHeat => 1.5); subject=Susceptibility{(:x, :x)}
    )
end

@testset "§8a load guard: >1 family slot rejected (needs §8b unification)" begin
    # two family-generic slots need cross-slot index unification (§8b) — rejected at load
    @test_throws Exception @eval @relation :test _TwoFam(
        a::Susceptibility, b::Conductivity
    ) = a - b
    # a plain abstract (non-family, non-concrete) slot is still rejected
    @test_throws Exception @eval @relation :test _AbstractSlotFD(x::AbstractQuantity) = x
end
