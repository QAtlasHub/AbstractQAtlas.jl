# @bound — the three declaration forms, the roles they record, and the guards.
#
# The load-bearing test here is the FIRST one: every declared bound must return a
# NUMBER as its slack.  During this migration a statement-form declaration whose
# body reached the macro wrapped in a `:block` fell through to the bare-slack
# branch, so the kernel evaluated `τ >= π/(2ΔE)` instead of measuring it and
# returned `true` — and `true ≥ -atol` holds, so the bound passed on every input,
# violations included.
#
# The per-relation tests DID catch it (measured: reintroducing the bug fails 15
# assertions across test_entanglement.jl and test_inequalities.jl), because every
# bound in the suite today has a saturation or fabricated-violation case pinning
# its slack to a number.  The sweep below is the STRUCTURAL version of that: it
# holds for a bound declared tomorrow with no such case, where the per-relation
# net has no thread.  A bound whose only failure mode is silence needs a check
# that does not depend on someone remembering to write one.

using AbstractQAtlas
using Test
using AbstractQAtlas:
    AbstractInequality,
    all_relations,
    bound_direction,
    bounded_slot,
    bounding_constant,
    bounding_slot,
    bounds_on,
    check,
    quantities,
    residual,
    slack,
    solve,
    variable_slots,
    variables

_bounds() = filter(r -> r isa AbstractInequality, all_relations())

# `@bound` pushes into the global registry as a load-time side effect.  The probes
# below are declarations under test, not physics, so each is withdrawn right after
# it is declared — otherwise `test_interface.jl`'s relation count (an invariant
# worth keeping exact) would move every time this file grows a case.  The struct
# and its methods survive; only registry membership is dropped.
function _unregister!(::Type{T}) where {T}
    return filter!(r -> !(r isa T), AbstractQAtlas._RELATION_REGISTRY)
end

# One probe per FORM, so each expansion path is exercised here and not only at
# package load: the physics declarations all expand during precompilation, where
# a form that stopped working would take the whole package down rather than fail
# a test — which is louder, but tells you nothing about which form broke.
@bound :test _BoolSlackProbe(a, b) = a - b                            # bare slack
@bound :test _BareSlackBound(x, y) = x - y                            # bare slack
@bound :test _CmpProbe(u <= u_max::SpecificHeat)                      # comparison, typed bounding
@bound :test _CmpConstProbe(w >= 2)                                   # comparison, constant bounding
@bound :test _StmtProbe(p, q) = p >= q / 2                            # statement, expression bounding
foreach(
    _unregister!, (_BoolSlackProbe, _BareSlackBound, _CmpProbe, _CmpConstProbe, _StmtProbe)
)

@testset "every declared bound measures its slack — never evaluates it" begin
    bs = _bounds()
    @test length(bs) >= 20            # not vacuous: the registry really is populated
    for r in bs
        # a positive, finite probe: no relation here divides by or logs a zero at 1
        kw = NamedTuple(v => 1.0 for v in variables(r))
        s = residual(r; kw...)
        @test !(s isa Bool)
        @test s isa Number
    end
end

@testset "a Bool slack is rejected loudly, not silently always-passed" begin
    # The guard, exercised directly: a kernel that returns a Bool would make
    # `slack ≥ -atol` true for every input.
    @test residual(_BoolSlackProbe(); a=2.0, b=1.0) == 1.0
    # forcibly install a comparing kernel and confirm the verbs refuse it
    AbstractQAtlas._residual(::_BoolSlackProbe; a, b, kw...) = a >= b
    @test_throws ErrorException residual(_BoolSlackProbe(); a=2.0, b=1.0)
    @test_throws ErrorException check(_BoolSlackProbe(); a=0.0, b=1.0)
end

@testset "comparison form records subject, bound and direction" begin
    # `bounded OP bounding`, both slots, one identity-typed
    @test bounded_slot(LiebRobinsonBound()) === :v
    @test bounding_slot(LiebRobinsonBound()) === :v_LR
    @test bounding_constant(LiebRobinsonBound()) === nothing
    @test bound_direction(LiebRobinsonBound()) === :upper
    @test variables(LiebRobinsonBound()) == (:v, :v_LR)
    # the bounded side stays UNTYPED while the bounding side keys a quantity —
    # the shape that lets a bound be stated before the bounded observable has a type
    @test variable_slots(LiebRobinsonBound()) ==
        ((:v, nothing), (:v_LR, LiebRobinsonVelocity))

    # `bounded OP constant`
    @test bounded_slot(SpecificHeatPositivity()) === :Cv
    @test bounding_slot(SpecificHeatPositivity()) === nothing
    @test bounding_constant(SpecificHeatPositivity()) == 0
    @test bound_direction(SpecificHeatPositivity()) === :lower
    @test bounding_constant(IoffeRegel()) == 1
    @test bound_direction(IoffeRegel()) === :lower
end

@testset "statement form records what it honestly can" begin
    # a bare-variable subject with an EXPRESSION bounding it: subject recorded,
    # bounding slot `nothing` (there is no single one) — not guessed
    @test bounded_slot(Subadditivity()) === :S_AB
    @test bounding_slot(Subadditivity()) === nothing
    @test bounding_constant(Subadditivity()) === nothing
    @test bound_direction(Subadditivity()) === :upper

    @test bounded_slot(MandelstamTammBound()) === :τ
    @test bound_direction(MandelstamTammBound()) === :lower

    # neither side is a single variable — no subject to claim
    @test bounded_slot(StrongSubadditivity()) === nothing
    @test bounding_slot(StrongSubadditivity()) === nothing
    @test bound_direction(StrongSubadditivity()) === :upper
end

@testset "bare-slack form: no direction is claimed" begin
    @test _BareSlackBound() isa AbstractInequality
    @test bound_direction(_BareSlackBound()) === :slack
    @test bounded_slot(_BareSlackBound()) === nothing
    @test bounding_slot(_BareSlackBound()) === nothing
    @test slack(_BareSlackBound(); x=3.0, y=1.0) == 2.0
    @test check(_BareSlackBound(); x=3.0, y=1.0)
    @test !check(_BareSlackBound(); x=1.0, y=3.0)
end

@testset "each form expands here, not only at package load" begin
    # comparison form with a TYPED bounding slot: the slot list is rebuilt from the
    # comparison, left to right, and the type reaches `variable_slots`
    @test variables(_CmpProbe()) == (:u, :u_max)
    @test variable_slots(_CmpProbe()) == ((:u, nothing), (:u_max, SpecificHeat))
    @test bounded_slot(_CmpProbe()) === :u
    @test bounding_slot(_CmpProbe()) === :u_max
    @test bound_direction(_CmpProbe()) === :upper
    @test slack(_CmpProbe(); u=1.0, u_max=3.0) == 2.0
    @test quantities(_CmpProbe()) == (SpecificHeat,)

    # comparison against a non-zero constant: the constant is recorded AND used
    @test bounding_constant(_CmpConstProbe()) == 2
    @test variables(_CmpConstProbe()) == (:w,)
    @test slack(_CmpConstProbe(); w=5) === 3
    @test !check(_CmpConstProbe(); w=1)

    # statement form with an expression on the bounding side
    @test bounded_slot(_StmtProbe()) === :p
    @test bounding_slot(_StmtProbe()) === nothing
    @test bound_direction(_StmtProbe()) === :lower
    @test slack(_StmtProbe(); p=3.0, q=4.0) == 1.0
end

@testset "the role validator rejects a slot that names nothing" begin
    # Defensive: the macro cannot produce this, but a hand-written `bounded_slot`
    # method can, and a role pointing at no variable makes every consumer reading
    # it silently wrong.  Exercised by calling the validator directly.
    AbstractQAtlas.bounded_slot(::_BareSlackBound) = :not_a_variable
    @test_throws ErrorException AbstractQAtlas._validate_bound_roles(_BareSlackBound())
    AbstractQAtlas.bounded_slot(::_BareSlackBound) = nothing          # restore
    @test AbstractQAtlas._validate_bound_roles(_BareSlackBound()) === nothing
end

@testset "an equality relation claims no bound roles" begin
    @test bound_direction(Rushbrooke()) === nothing
    @test bounded_slot(Rushbrooke()) === nothing
    @test bounding_slot(Rushbrooke()) === nothing
    @test bounding_constant(Rushbrooke()) === nothing
end

@testset "declaration guards" begin
    # strict operators would be DECLARED strict and CHECKED non-strictly
    @test_throws Exception @eval @bound :test _Strict(a < b)
    @test_throws Exception @eval @bound :test _StrictStmt(a, b) = a < b
    # the subject goes on the left; a constant there is a mis-stated bound
    @test_throws Exception @eval @bound :test _ConstLeft(0 <= x)
    # an expression side in the comparison form leaves its variables undeclared
    @test_throws Exception @eval @bound :test _ExprSide(a <= b + c)
    # one statement per bound
    @test_throws Exception @eval @bound :test _Two(a <= b, c <= d)
end

@testset "bounds_on answers 'what bounds this?' by ROLE, not membership" begin
    @test SpecificHeatPositivity() in bounds_on(SpecificHeat)
    # family-declared, so a concrete component finds it (§8a)
    @test SusceptibilityPositivity() in bounds_on(Susceptibility{(:z, :z)})
    # The sharp case: `LiebRobinsonVelocity` appears in BOTH queries, but for
    # different relations, because its ROLE differs in each.
    #   - in `LiebRobinsonBound(v <= v_LR)` it DOES the bounding, so the role-aware
    #     query must not return that relation...
    #   - ...while the role-blind `relations_constraining` does;
    #   - and in `VelocityPositivity(v >= 0)` it IS the subject, so `bounds_on`
    #     returns that one.
    # A `bounds_on` that ignored roles would return both; one that ignored the
    # `EachOf` quantifier would return neither.
    lrv_bounds = bounds_on(LiebRobinsonVelocity)
    @test !(LiebRobinsonBound() in lrv_bounds)
    @test VelocityPositivity() in lrv_bounds
    @test LiebRobinsonBound() in AbstractQAtlas.relations_constraining(LiebRobinsonVelocity)
end

@testset "slacks are unchanged by the migration" begin
    # Each expected value is the PRE-migration slack expression, written out here
    # independently — so a re-derivation of the same formula from the new
    # declaration cannot agree with itself.
    S_A, S_B, S_AB, S_BC, S_ABC = 0.6, 0.7, 1.1, 1.2, 1.4
    @test slack(Subadditivity(); S_A, S_B, S_AB) ≈ S_A + S_B - S_AB
    @test slack(StrongSubadditivity(); S_AB, S_BC, S_ABC, S_B) ≈ S_AB + S_BC - S_ABC - S_B
    @test slack(WeakMonotonicity(); S_AB, S_BC, S_A, S_C=0.5) ≈ S_AB + S_BC - S_A - 0.5
    @test slack(ArakiLieb(); S_AB, S_A, S_B) ≈ S_AB - abs(S_A - S_B)
    @test slack(MaxEntropyBound(); S=0.4, log_d=log(2)) ≈ log(2) - 0.4
    @test slack(RenyiMonotonicity(); S_low=0.9, S_high=0.4) ≈ 0.9 - 0.4
    @test slack(EntropyMixingConcavity(); S_mix=0.9, S_avg=0.4) ≈ 0.9 - 0.4
    @test slack(HolevoMixingBound(); S_avg=0.4, H_weights=log(2), S_mix=0.9) ≈
        0.4 + log(2) - 0.9
    @test slack(MeasurementEntropyIncrease(); S_meas=0.9, S=0.4) ≈ 0.9 - 0.4
    @test slack(Monogamy(); τ_ABC=1.0, τ_AB=0.3, τ_AC=0.2) ≈ 1.0 - 0.3 - 0.2
    @test slack(RobertsonUncertainty(); ΔA=0.8, ΔB=0.9, comm=1.0) ≈ 0.8 * 0.9 - 1.0 / 2
    @test slack(MandelstamTammBound(); τ=2π, ΔE=0.5) ≈ 2π - π / (2 * 0.5)
    @test slack(MargolusLevitinBound(); τ=2π, E_above=0.5) ≈ 2π - π / (2 * 0.5)
    @test slack(LiebRobinsonBound(); v=2.0, v_LR=3.0) ≈ 3.0 - 2.0
    @test slack(CTheorem(); c_UV=1.0, c_IR=0.5) ≈ 1.0 - 0.5
    @test slack(JarzynskiSecondLaw(); W_avg=1.0, ΔF=0.6) ≈ 1.0 - 0.6
    @test slack(IoffeRegel(); kFℓ=3.0) ≈ 3.0 - 1
    @test slack(SpecificHeatPositivity(); Cv=0.7) ≈ 0.7
    @test slack(CompressibilityPositivity(); κT=0.7) ≈ 0.7
    @test slack(SusceptibilityPositivity(); χT=0.7) ≈ 0.7
    @test slack(EntropyNonNegativity(); S=0.7) ≈ 0.7
    @test slack(RelativeEntropyNonNegativity(); S_rel=0.7) ≈ 0.7
end

@testset "a constant-bounded slack keeps the exact-arithmetic contract" begin
    # `Cv >= 0` must generate the bare subject, not `Cv - 0`: a Rational in stays
    # a Rational out, and an Int slack stays an Int.
    r = slack(SpecificHeatPositivity(); Cv=3//2)
    @test r === 3//2
    @test slack(EntropyNonNegativity(); S=0//1) === 0//1
    @test slack(IoffeRegel(); kFℓ=3) === 2
end

@testset "declared direction agrees with the criterion it generates" begin
    # An :upper bound must FAIL when its subject exceeds the bounding value, and a
    # :lower bound when the subject falls below it — the direction metadata and the
    # generated slack cannot drift apart without this failing.
    for r in _bounds()
        d = bound_direction(r)
        bd, bg = bounded_slot(r), bounding_slot(r)
        (bd === nothing || bg === nothing) && continue
        base = NamedTuple(v => 1.0 for v in variables(r))
        violating = d === :upper ? merge(base, (; bd => 2.0)) : merge(base, (; bd => 0.5))
        satisfying = d === :upper ? merge(base, (; bd => 0.5)) : merge(base, (; bd => 2.0))
        @test !check(r; violating...)
        @test check(r; satisfying...)
    end
end
