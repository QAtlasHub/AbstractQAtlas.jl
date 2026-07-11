# Fourier / conjugate-representation structure: real ↔ momentum, time ↔
# frequency, and the quantity pairs related by the transform.

using AbstractQAtlas
using AbstractQAtlas:
    representation, fourier_conjugate, fourier_conjugate_quantity, fourier_pair

@testset "conjugate representations form an involution" begin
    @test fourier_conjugate(RealSpace()) === MomentumSpace()
    @test fourier_conjugate(MomentumSpace()) === RealSpace()
    @test fourier_conjugate(TimeDomain()) === FrequencyDomain()
    @test fourier_conjugate(FrequencyDomain()) === TimeDomain()
    for r in (RealSpace(), MomentumSpace(), TimeDomain(), FrequencyDomain())
        @test fourier_conjugate(fourier_conjugate(r)) === r
    end
end

@testset "each quantity's representation" begin
    @test representation(SpinCorrelation(:z, :z)) == (RealSpace(),)
    @test representation(StaticStructureFactor()) == (MomentumSpace(),)
    @test representation(DynamicalCorrelation(:x, :y)) == (RealSpace(), TimeDomain())
    @test representation(DynamicalStructureFactor()) == (MomentumSpace(), FrequencyDomain())
    @test representation(DynamicalSusceptibility(:x, :y)) ==
        (MomentumSpace(), FrequencyDomain())
    @test representation(RetardedGreensFunction()) == (MomentumSpace(), FrequencyDomain())
    # AC conductivity in (q, ω), its current–current correlation in (r, t)
    @test representation(DynamicalConductivity(:x, :y)) ==
        (MomentumSpace(), FrequencyDomain())
    @test representation(CurrentCorrelation(:x, :y)) == (RealSpace(), TimeDomain())
    @test representation(CurrentNoise(:x, :y)) == (MomentumSpace(), FrequencyDomain())
    # global thermodynamic quantities carry no space/time representation
    @test representation(Energy()) == ()
    @test representation(FreeEnergy()) == ()
    @test representation(SpecificHeat()) == ()
end

@testset "Fourier-conjugate quantity + pairing" begin
    @test fourier_conjugate_quantity(StaticStructureFactor()) === SpinCorrelation
    @test fourier_conjugate_quantity(DynamicalStructureFactor()) === DynamicalCorrelation
    @test fourier_conjugate_quantity(DynamicalCorrelation(:x, :y)) ===
        DynamicalStructureFactor
    @test fourier_conjugate_quantity(Energy()) === nothing

    # current channel mirrors the spin channel: S^j(q,ω) ↔ ⟨jj⟩(r,t)
    @test fourier_conjugate_quantity(CurrentNoise(:x, :y)) === CurrentCorrelation
    @test fourier_conjugate_quantity(CurrentCorrelation(:x, :y)) === CurrentNoise
    @test fourier_pair(CurrentNoise(:x, :y), CurrentCorrelation(:x, :y))

    # spatial FT: S(q) ↔ ⟨S S⟩(r); space-time FT: S(q,ω) ↔ ⟨A A⟩(r,t)
    @test fourier_pair(StaticStructureFactor(), SpinCorrelation(:z, :z))
    @test fourier_pair(DynamicalStructureFactor(), DynamicalCorrelation(:x, :y))
    # not conjugates: different number of representations, or same space
    @test !fourier_pair(StaticStructureFactor(), DynamicalStructureFactor())
    @test !fourier_pair(StaticStructureFactor(), StaticStructureFactor())
    @test !fourier_pair(Energy(), Energy())
    # the representations of a conjugate pair are elementwise conjugate
    rS = representation(DynamicalStructureFactor())
    rG = representation(DynamicalCorrelation(:x, :y))
    @test all(fourier_conjugate(x) === y for (x, y) in zip(rS, rG))
end
