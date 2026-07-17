# The declarative interface itself: registry, introspection, generic
# affine solve (incl. its refusal on non-affine variables), and the
# one-call adoption surface.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

# a deliberately non-affine relation (kernel wired by hand, NOT registered,
# so registry counts stay untouched) to exercise the generic-solve guard
struct _NonAffineDemo <: AbstractRelation end
AbstractQAtlas._residual(::_NonAffineDemo; a, b, _extra...) = a - b^2
AbstractQAtlas.variables(::_NonAffineDemo) = (:a, :b)
AbstractQAtlas.domain(::_NonAffineDemo) = :test_only

@testset "registry + traits" begin
    rels = all_relations()
    @test length(rels) == 105         # universal-only: model-specific (spin glass, Drude mobility, single-band Hall) moved to QAtlas
    @test allunique(typeof.(rels))
    @test length(all_relations(; domain=:scaling)) == 5
    @test length(all_relations(; domain=:thermodynamic)) == 15
    @test length(all_relations(; domain=:fundamental)) == 6
    @test length(all_relations(; domain=:topology)) == 3
    @test length(all_relations(; domain=:spectral)) == 10
    @test length(all_relations(; domain=:keldysh)) == 6
    @test length(all_relations(; domain=:transport)) == 18
    @test length(all_relations(; domain=:quantum)) == 9
    @test length(all_relations(; domain=:ensemble)) == 2
    @test length(all_relations(; domain=:entanglement)) == 25
    @test length(all_relations(; domain=:wick)) == 2
    @test length(all_relations(; domain=:cft)) == 4
    @test isempty(all_relations(; domain=:spinglass))   # model-specific — lives in QAtlas now
    @test variables(Widom()) == (:β, :γ, :δ)
    @test variables(SpecificHeatFDT()) == (:C, :var_E, :β)   # N optional, not listed
    @test domain(Rushbrooke()) == :scaling
    @test domain(TKNN()) == :topology
end

@testset "generic affine solve: exactness + refusal" begin
    # generic solve is exact for Rational data (no hand-written rearrangement)
    @test solve(Widom(), Val(:δ); β=1//8, γ=7//4) == 15//1
    @test solve(FreeEnergyLegendre(), Val(:S); F=1//2, U=3//2, β=2//1) == 2//1
    @test solve(TKNN(), Val(:σxy); C=2) == 2
    # non-affine variable is REFUSED, not silently mis-solved:
    err = try
        solve(_NonAffineDemo(), Val(:b); a=4.0)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("not affine", err.msg)
    # …while its affine variable solves generically (a = b²)
    @test solve(_NonAffineDemo(), Val(:a); b=3) == 9
    # solving FOR the inverse temperature: the β-or-T demand is waived
    # for the target itself, so SpecificHeatFDT (quadratic in β) hits the
    # affinity guard — while the β-linear susceptibility solves exactly:
    @test_throws ErrorException solve(SpecificHeatFDT(), Val(:β); C=1.0, var_E=1.0)
    @test solve(SusceptibilityFDT(), Val(:β); χ=1//2, var_M=2//1, N=4) == 1//1
    # a variable the relation does not depend on is refused too
    @test_throws ErrorException solve(Rushbrooke(), Val(:η); α=0//1, β=1//8, γ=7//4)
end

@testset "β-or-T normalization at every verb" begin
    @test residual(SusceptibilityFDT(); χ=1.0, var_M=2.0, T=2.0) ==
        residual(SusceptibilityFDT(); χ=1.0, var_M=2.0, β=0.5)
    @test_throws ErrorException residual(SusceptibilityFDT(); χ=1.0, var_M=2.0)
    @test_throws ErrorException residual(
        SusceptibilityFDT(); χ=1.0, var_M=2.0, β=0.5, T=2.0
    )
    # relations without β are untouched by a stray T… (not applicable) —
    # they simply never list it:
    @test :β ∉ variables(TKNN())
end

@testset "one-call adoption: applicable_relations / relation_report / check_all" begin
    # exponent table: exactly the four scaling laws apply
    ising = (α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4, d=2)
    app = applicable_relations(ising)
    @test length(app) == 4
    @test all(r -> domain(r) == :scaling, app)
    @test check_all(ising)                      # atol = 0: exact gate
    rep = relation_report(ising)
    @test all(row -> row.pass, rep)
    @test all(row -> row.residual isa Rational, rep)   # exactness through the sweep

    # measured thermodynamics: exactly SpecificHeatFDT applies, optional N honored
    β = 0.8
    varE = 1.7
    N = 4
    C = β^2 * varE / N
    thermo = (C=C, var_E=varE, β=β, N=N)
    app2 = applicable_relations(thermo)
    @test length(app2) == 1
    @test app2[1] isa SpecificHeatFDT
    @test check_all(thermo; atol=1e-14)
    # T-form data normalizes too
    @test check_all((C=C, var_E=varE, T=1 / β, N=N); atol=1e-14)

    # empty match is false, never a silent green
    @test !check_all((; unrelated=1.0))

    # a wrong value is caught by the sweep
    @test !check_all((; ising..., γ=3//2))
end

@testset "wrappers delegate to the registry sweep" begin
    ising = (α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4)
    @test exponents_consistent(ising; d=2)
    r = exponent_residuals(ising; d=2)
    @test keys(r) == (:rushbrooke, :widom, :fisher, :josephson)
    @test all(iszero, values(r))
end
