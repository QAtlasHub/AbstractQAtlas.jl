# Gap ↔ correlation length, and the dynamical exponent.
#
# ξ = v/Δ is verified against the real-space decay of a relativistic
# dispersion; the dynamical exponent is read off a synthetic Δ(ξ) power
# law.  This file also pins the generic-solve guard against a probe that
# blows up (ξ = v/Δ is non-affine in Δ).

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

@testset "CorrelationLengthGap: ξ = v/Δ" begin
    for (v, Δ) in ((1.0, 0.5), (2.3, 0.1), (0.7, 1.4))
        @test check(CorrelationLengthGap(); ξ=v / Δ, v=v, Δ=Δ, atol=1e-13)
        # affine in ξ and v — those solve; Δ (1/Δ) does not
        @test solve(CorrelationLengthGap(), Val(:ξ); v=v, Δ=Δ) ≈ v / Δ
        @test solve(CorrelationLengthGap(), Val(:v); ξ=v / Δ, Δ=Δ) ≈ v
    end
    @test !check(CorrelationLengthGap(); ξ=1.0, v=1.0, Δ=1.0 + 1e-3)
end

@testset "generic solve refuses a non-affine variable that blows up at a probe" begin
    # ξ = v/Δ is 1/Δ: solving for Δ would probe Δ = 0 (⇒ Inf); the guard
    # must REFUSE, not silently return NaN.
    err = try
        solve(CorrelationLengthGap(), Val(:Δ); ξ=2.0, v=1.0)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("not affine", err.msg)
end

@testset "DynamicalScaling: d(lnΔ)/d(lnξ) = −z reads off the dynamical exponent" begin
    for z in (1.0, 2.0, 0.5)
        A = 3.0
        Δ(ξ) = A * ξ^(-z)                      # synthetic gap–ξ power law
        ξ0, h = 10.0, 1e-5
        slope = (log(Δ(ξ0 * exp(h))) - log(Δ(ξ0 * exp(-h)))) / (2h)   # d lnΔ/d lnξ
        @test check(DynamicalScaling(); dlogΔ_dlogξ=slope, z=z, atol=1e-6)
        @test solve(DynamicalScaling(), Val(:z); dlogΔ_dlogξ=slope) ≈ z atol = 1e-6
    end
    # z = 1 is the Lorentz-invariant (relativistic) value
    @test check(DynamicalScaling(); dlogΔ_dlogξ=-1.0, z=1.0)
end

@testset "a NAMED velocity reaches ξ = v/Δ (the family slot is the point)" begin
    # Before `Velocity{K}`, `FermiVelocity` was its own struct, so `typeof` made it a
    # different bag key from `Velocity` and this relation could not see it — every
    # atlas hub that knows its Fermi or Luttinger velocity rather than an anonymous
    # "velocity" was cut out of ξ = v/Δ and of the two CFT finite-size relations.
    v, Δ = 2.5, 0.5                                   # exact in binary: ξ = 5.0
    b = bag(FermiVelocity() => v, MassGap => Δ, CorrelationLength => v / Δ)

    # The #802 metric, stated directly: a named velocity is now CONSTRAINED by the
    # three relations whose `v` slot is the family. It was constrained by none.
    constrained = Set(typeof(r) for r in relations_constraining(FermiVelocity))
    @test CorrelationLengthGap in constrained
    @test FiniteSizeGap in constrained
    @test CasimirCentralCharge in constrained

    @test any(r -> r isa CorrelationLengthGap, applicable_relations(b))
    rows = [r for r in relation_report(b) if r.relation isa CorrelationLengthGap]
    @test length(rows) == 1
    @test rows[1].subject.type === Velocity{:fermi}   # the SUBJECT is the Fermi velocity
    @test rows[1].pass

    # One velocity in the bag ⇒ nothing to choose between ⇒ no `subject` ceremony,
    # which is also what lets `solve` (which passes none) traverse the relation.
    @test check(CorrelationLengthGap(), b; atol=1e-12)
    @test solve(CorrelationLengthGap(), CorrelationLength, b) ≈ v / Δ
end

@testset "two velocity kinds in one bag is a real ambiguity, and refuses" begin
    b = bag(
        FermiVelocity() => 2.0,
        LuttingerVelocity() => 3.0,
        MassGap => 0.5,
        CorrelationLength => 4.0,
    )
    err = try
        check(CorrelationLengthGap(), b)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("pass `subject`", err.msg)          # never silently pick one
    @test check(CorrelationLengthGap(), b; subject=Velocity{:fermi}, atol=1e-12)
    # ...and the report enumerates BOTH kinds rather than collapsing them
    rows = [r for r in relation_report(b) if r.relation isa CorrelationLengthGap]
    @test Set(r.subject.type for r in rows) ==
        Set((Velocity{:fermi}, Velocity{:luttinger}))
end
