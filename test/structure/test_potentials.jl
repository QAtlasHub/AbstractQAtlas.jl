# test/structure/test_potentials.jl — the thermodynamic-potential structure and the
# Maxwell relations DERIVED from it.  The key assertion: each hand-written Maxwell
# @relation is a CONSEQUENCE of the potential structure (its residual equals the
# structure-derived one for any inputs), so the four are instances of ONE structural
# fact — the equality of the mixed second partials — not four independent axioms.

using AbstractQAtlas
using Test
const AQ = AbstractQAtlas

@testset "thermodynamic_potentials: the four differentials" begin
    pots = thermodynamic_potentials()
    @test length(pots) == 4
    @test Set(p.name for p in pots) == Set((:U, :F, :H, :G))
    F = only(p for p in pots if p.name === :F)
    @test (F.x.variable, F.x.conjugate, F.x.sign) == (:T, :S, -1)   # dF ⊃ −S dT
    @test (F.y.variable, F.y.conjugate, F.y.sign) == (:V, :p, -1)   # dF ⊃ −p dV
end

@testset "maxwell_relation: derived from the potential structure" begin
    byname = Dict(p.name => maxwell_relation(p) for p in thermodynamic_potentials())
    # F(T,V): ∂S/∂V =  ∂p/∂T   (coeff +1)
    @test byname[:F].lhs == (:S, :V) && byname[:F].rhs == (:p, :T) && byname[:F].coeff == 1
    # G(T,p): ∂S/∂p = −∂V/∂T   (coeff −1)
    @test byname[:G].lhs == (:S, :p) && byname[:G].rhs == (:V, :T) && byname[:G].coeff == -1
    # U(S,V): ∂T/∂V = −∂p/∂S   (coeff −1)
    @test byname[:U].lhs == (:T, :V) && byname[:U].rhs == (:p, :S) && byname[:U].coeff == -1
    # H(S,p): ∂T/∂p =  ∂V/∂S   (coeff +1)
    @test byname[:H].lhs == (:T, :p) && byname[:H].rhs == (:V, :S) && byname[:H].coeff == 1
    # zero residual when the mixed partials commute; exact (Rational in ⇒ Rational out)
    r = maxwell_residual(byname[:F]; derivs=Dict((:S, :V) => 2 // 1, (:p, :T) => 2 // 1))
    @test r == 0 // 1 && r isa Rational
    # a missing derivative errors loudly (both slots checked)
    @test_throws ErrorException maxwell_residual(byname[:F]; derivs=Dict((:S, :V) => 1.0))
    @test_throws ErrorException maxwell_residual(byname[:F]; derivs=Dict((:p, :T) => 1.0))
end

@testset "the four Maxwell @relations are CONSEQUENCES of the potential structure" begin
    # each hand-written Maxwell relation + its variables in (∂cₓ/∂y, ∂c_y/∂x) order
    hand = Dict(
        :F => (MaxwellHelmholtz(), (:dS_dV, :dp_dT)),
        :G => (MaxwellGibbs(), (:dS_dp, :dV_dT)),
        :U => (MaxwellInternal(), (:dT_dV, :dp_dS)),
        :H => (MaxwellEnthalpy(), (:dT_dp, :dV_dS)),
    )
    for p in thermodynamic_potentials()
        m = maxwell_relation(p)
        rel, (v1, v2) = hand[p.name]
        # for ANY derivative values, the structure-derived residual EQUALS the
        # hand-written relation's residual — they are the same identity (the commuting
        # mixed partials), so the four @relations are its consequences, not axioms.
        for (a, b) in ((3 // 5, 7 // 5), (-1 // 2, 2 // 1), (1.0, 1.0), (0.3, -0.8))
            r_struct = maxwell_residual(m; derivs=Dict(m.lhs => a, m.rhs => b))
            r_hand = residual(rel; NamedTuple{(v1, v2)}((a, b))...)
            @test r_struct == r_hand
        end
    end
end

@testset "show covers all branches" begin
    pots = thermodynamic_potentials()
    F = only(p for p in pots if p.name === :F)   # both terms sign<0  (− branch of trm)
    U = only(p for p in pots if p.name === :U)   # has a +T dS term   (+ branch of trm)
    G = only(p for p in pots if p.name === :G)   # coeff −1           (coeff<0 branch)
    H = only(p for p in pots if p.name === :H)   # coeff +1           (coeff>0 branch)
    @test occursin("F(T,V)", sprint(show, F)) && occursin("U(S,V)", sprint(show, U))
    @test occursin("MaxwellRelation(G", sprint(show, maxwell_relation(G)))
    @test occursin("MaxwellRelation(H", sprint(show, maxwell_relation(H)))
end
