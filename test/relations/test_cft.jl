# CFT finite-size relations: read the central charge and scaling
# dimensions off synthetic finite-size spectra, the way an MPS/ED
# calculation would.  Independent check = the extracted universal number
# is the SAME across system sizes (the hallmark of a universal amplitude).

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

@testset "CasimirCentralCharge: e₀(L) = e_∞ − πcv/(6L²) reads off c" begin
    c, v = 1 / 2, 1.0                      # 2D Ising universality
    e∞ = -1.2732
    e0(L) = e∞ - π * c * v / (6 * L^2)     # synthetic finite-size energy density
    for L in (16, 32, 64, 128)
        dE = e0(L) - e∞                    # the finite-size correction
        @test check(CasimirCentralCharge(); dE=dE, c=c, v=v, L=L, atol=1e-13)
        # c recovered identically at every size (universal amplitude)
        @test solve(CasimirCentralCharge(), Val(:c); dE=dE, v=v, L=L) ≈ c
    end
    # a wrong c fails
    L = 32
    @test !check(CasimirCentralCharge(); dE=e0(L) - e∞, c=1.0, v=v, L=L)
end

@testset "FiniteSizeGap: E_x(L) − E₀(L) = 2πvx/L reads off the scaling dimension" begin
    v = 1.3
    for x in (1 / 8, 1.0, 15 / 8)          # Ising σ, ε, and a descendant
        gaps = Dict(L => 2π * v * x / L for L in (16, 32, 64))
        for (L, gap) in gaps
            @test check(FiniteSizeGap(); gap=gap, x=x, v=v, L=L, atol=1e-12)
            @test solve(FiniteSizeGap(), Val(:x); gap=gap, v=v, L=L) ≈ x
        end
        # the extracted dimension is L-independent (universal)
        xs = [solve(FiniteSizeGap(), Val(:x); gap=gaps[L], v=v, L=L) for L in (16, 32, 64)]
        @test all(≈(x), xs)
    end
end

@testset "one-call sweep: finite-size data → central charge + dimension" begin
    # a bag of measured finite-size numbers picks up exactly the cft relations
    data = (
        dE=(-π * 0.5 * 1.0 / (6 * 64^2)),
        c=0.5,
        v=1.0,
        L=64,
        gap=2π * 1.0 * 0.125 / 64,
        x=0.125,
    )
    rep = relation_report(data; atol=1e-12, domain=:cft)
    @test length(rep) == 2
    @test all(row -> row.pass, rep)
end

@testset "Cardy density of states + Zamolodchikov c-theorem" begin
    using AbstractQAtlas: check, solve, slack, AbstractInequality
    # ln ρ(Δ) = 2π√(cΔ/6): Ising c=1/2 at Δ=3 ⇒ ln ρ = 2π√(1/4) = π
    @test check(
        CardyDensityOfStates(); ln_ρ=2π * sqrt((1 / 2) * 3 / 6), c=1 / 2, Δ=3, atol=1e-12
    )
    @test solve(CardyDensityOfStates(), Val(:ln_ρ); c=1 / 2, Δ=3) ≈ 2π * sqrt(1 / 4)
    @test solve(CardyDensityOfStates(), Val(:ln_ρ); c=1 / 2, Δ=3) ≈ π atol = 1e-12
    # c-theorem c_UV ≥ c_IR: tricritical→critical Ising 7/10 → 1/2 (an allowed flow)
    @test CTheorem() isa AbstractInequality
    @test check(CTheorem(); c_UV=7 // 10, c_IR=1 // 2)
    @test slack(CTheorem(); c_UV=1 // 2, c_IR=1 // 2) == 0 // 1     # a fixed point (no flow)
    @test !check(CTheorem(); c_UV=1 // 2, c_IR=7 // 10)            # c increasing ⇒ forbidden
end

@testset "type-keyed: cft" begin
    @test quantities(CasimirCentralCharge()) == (CentralCharge,)
    @test Set(quantities(FiniteSizeGap())) == Set((MassGap, ScalingDimension))
    @test Set(quantities(CardyDensityOfStates())) == Set((CentralCharge, ScalingDimension))
    # finite-size gap Δ = 2π v x / L via bag (MassGap, ScalingDimension typed; v, L supplied)
    v, L, x = 1.0, 10.0, 0.125
    @test check(
        FiniteSizeGap(),
        bag(MassGap => 2π * v * x / L, ScalingDimension => x);
        v=v,
        L=L,
        atol=1e-12,
    )
    # CTheorem c_UV ≥ c_IR stays symbol-keyed — c_UV/c_IR are two CentralCharge
    # instances (UV/IR), the multi-instance-of-one-type case (Phase-2 decoration)
    @test isempty(variable_types(CTheorem()))
end
