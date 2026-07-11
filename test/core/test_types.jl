# Core type vocabulary: subtyping, BC size semantics, fetch fallback.

using AbstractQAtlas
using AbstractQAtlas: _bc_size, fetch

@testset "type tree" begin
    @test Universality(:Ising) isa AbstractQAtlasModel
    @test Universality(:Ising) === Universality{:Ising}()
    @test Infinite() isa BoundaryCondition
    @test OBC(4) isa BoundaryCondition
    @test PBC(4) isa BoundaryCondition
    @test CriticalExponents() isa AbstractQuantity
    @test SpecificHeat() isa AbstractThermalPotential
    @test SpecificHeat() isa AbstractQuantity
    @test Magnetization(:z) isa AbstractMagnetization
    @test Susceptibility(:z, :z) isa AbstractSusceptibility
    @test PartitionFunction() isa AbstractThermalPotential
    @test SpontaneousMagnetization() isa AbstractMagnetization
    @test CriticalTemperature() isa AbstractQuantity
    @test TopologicalInvariant() isa AbstractQuantity
end

@testset "tensor traits (indices / rank / spaces)" begin
    # honest index tuples, one entry per tensor slot
    @test indices(Magnetization(:x)) == (:x,)
    @test indices(Magnetization(:z)) == (:z,)
    @test indices(Susceptibility(:x, :x)) == (:x, :x)
    @test indices(Susceptibility(:z, :z)) == (:z, :z)
    @test indices(Susceptibility(:x, :y)) == (:x, :y)          # off-diagonal expressible
    # scalars: empty index tuple, rank 0
    @test indices(SpecificHeat()) == ()
    @test indices(CriticalExponents()) == ()
    @test tensor_rank(SpecificHeat()) == 0
    @test tensor_rank(Magnetization(:x)) == 1
    @test tensor_rank(Susceptibility(:x, :y)) == 2
    @test index_spaces(Susceptibility(:x, :y)) == (SpinAxis(), SpinAxis())
    @test index_spaces(RetardedGreensFunction()) == (OrbitalIndex(), OrbitalIndex())
    @test index_spaces(Conductivity(:x, :y)) == (SpatialDirection(), SpatialDirection())
end

@testset "nonlinear response order (higher-order susceptibility)" begin
    # χ⁽ⁿ⁾_{α;β₁…βₙ} = ∂ⁿM_α/∂hⁿ : n = length(indices) − 1
    @test response_order(Susceptibility(:x, :y)) == 1            # linear
    @test response_order(Susceptibility(:x, :y, :z)) == 2        # 2nd-order nonlinear
    @test response_order(Susceptibility(:x, :x, :x, :x)) == 3    # 3rd-order
    @test tensor_rank(Susceptibility(:x, :y, :z)) == 3           # rank = n + 1
    @test index_spaces(Susceptibility(:x, :y, :z)) == (SpinAxis(), SpinAxis(), SpinAxis())
    @test response_order(Conductivity(:x, :y, :z)) == 2          # nonlinear conductivity
    @test response_order(SpecificHeat()) == 0                    # not a response function
    # a susceptibility needs ≥2 indices (1 response + ≥1 field)
    @test_throws ErrorException Susceptibility(:x)
    @test_throws ErrorException Susceptibility{(:x, 2)}()        # non-symbol index
end

@testset "multi-time: frequency_arguments (dynamical nonlinear response)" begin
    # static χ⁽ⁿ⁾ = ∂ⁿM/∂hⁿ is the zero-frequency limit — no frequency args
    @test frequency_arguments(Susceptibility(:x, :y)) == 0
    @test frequency_arguments(Susceptibility(:x, :y, :z)) == 0
    # dynamical response: n-th order ⇒ n independent frequencies = multi-time
    @test frequency_arguments(DynamicalSusceptibility(:x, :y)) == 1       # χ(ω)
    @test frequency_arguments(DynamicalSusceptibility(:x, :y, :z)) == 2   # χ⁽²⁾(ω₁,ω₂)
    @test frequency_arguments(DynamicalSusceptibility(:x, :x, :x, :x)) == 3
    # the dynamical susceptibility carries the same tensor structure too
    @test response_order(DynamicalSusceptibility(:x, :y, :z)) == 2
    @test indices(DynamicalSusceptibility(:x, :y, :z)) == (:x, :y, :z)
    # one-frequency dynamical quantities
    @test frequency_arguments(RetardedGreensFunction()) == 1
    @test frequency_arguments(DynamicalStructureFactor()) == 1
    @test frequency_arguments(DensityOfStates()) == 1
    # static / instantaneous quantities: none
    @test frequency_arguments(Energy()) == 0
    @test frequency_arguments(SpecificHeat()) == 0
    @test frequency_arguments(Magnetization(:z)) == 0
    # Conductivity is the DC (static) response — the current-channel analogue
    # of the static Susceptibility — so it has NO frequency arguments at any
    # order.  Its AC counterpart DynamicalConductivity is frequency-resolved.
    @test frequency_arguments(Conductivity(:x, :y)) == 0
    @test frequency_arguments(Conductivity(:x, :y, :z)) == 0
    # AC conductivity σ⁽ⁿ⁾(ω₁…ωₙ): n-th order ⇒ n frequencies (like the
    # dynamical susceptibility), rank n+1 in SpatialDirection space
    @test frequency_arguments(DynamicalConductivity(:x, :y)) == 1
    @test frequency_arguments(DynamicalConductivity(:x, :y, :z)) == 2
    @test response_order(DynamicalConductivity(:x, :y, :z)) == 2
    @test index_spaces(DynamicalConductivity(:x, :y, :z)) ==
        (SpatialDirection(), SpatialDirection(), SpatialDirection())
    # its Kubo kernel, the current–current correlation, is n-time to match
    @test frequency_arguments(CurrentCorrelation(:x, :y)) == 1
    @test frequency_arguments(CurrentCorrelation(:x, :y, :z)) == 2
    @test index_spaces(CurrentCorrelation(:x, :y)) ==
        (SpatialDirection(), SpatialDirection())
    # both need ≥2 indices
    @test_throws ErrorException DynamicalConductivity(:x)
    @test_throws ErrorException CurrentCorrelation(:x)
end

@testset "Energy granularity" begin
    @test Energy() === Energy{:natural}()
    @test Energy(:total) === Energy{:total}()
    @test Energy(:per_site) === Energy{:per_site}()
    @test_throws ErrorException Energy(:bogus)
end

@testset "BC size semantics (N = 0 sentinel)" begin
    @test _bc_size(OBC(12), NamedTuple()) == 12
    @test _bc_size(PBC(8), NamedTuple()) == 8
    # sentinel: N from kwargs
    @test _bc_size(OBC(), (N=24,)) == 24
    @test _bc_size(PBC(; N=0), (N=6,)) == 6
    # bc.N wins over kwargs
    @test _bc_size(OBC(12), (N=99,)) == 12
    # unresolvable
    @test_throws ErrorException _bc_size(OBC(), NamedTuple())
    @test_throws ErrorException _bc_size(Infinite(), (N=4,))
end

@testset "fetch fallback throws informatively" begin
    struct _TTModel <: AbstractQAtlasModel end
    err = try
        fetch(_TTModel(), SpecificHeat(), Infinite())
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("_TTModel", err.msg)
    @test occursin("SpecificHeat", err.msg)
    @test occursin("Infinite", err.msg)
end
