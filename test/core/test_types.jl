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

@testset "transport family tensor structure" begin
    # currents are rank-1 vectors in SpatialDirection space
    @test tensor_rank(ElectricCurrent()) == 1
    @test tensor_rank(HeatCurrent()) == 1
    @test index_spaces(HeatCurrent()) == (SpatialDirection(),)
    # the transport coefficients are rank-2 tensors
    for Q in (DrudeWeight, ThermalConductivity, Thermopower, PeltierCoefficient)
        @test tensor_rank(Q(:x, :y)) == 2
        @test index_spaces(Q(:x, :y)) == (SpatialDirection(), SpatialDirection())
        @test indices(Q(:x, :y)) == (:x, :y)
        @test_throws ErrorException Q(:x)          # rank-2: needs exactly 2 indices
        @test_throws ErrorException Q(:x, :y, :z)
    end
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

@testset "TsallisEntropy carries its order too, and keys the same way" begin
    # Same shape as RenyiEntropy: `S_q` for two different `q` are two quantities, and a
    # bag keyed by type alone would keep only the last written. The relation
    # `TsallisEntropyMoment(Sq, moment, q)` still takes the order as a variable — the
    # field is what makes the BAG able to hold two, not what feeds the formula.
    @test TsallisEntropy(2).q == 2.0
    @test TsallisEntropy(2) != TsallisEntropy(3)
    @test_throws ArgumentError TsallisEntropy(0)
    @test_throws ArgumentError TsallisEntropy(-1)
    @test_throws ArgumentError TsallisEntropy(1)     # the von Neumann limit
    @test variable_support(TsallisEntropy(2)) == OrderSupport(2.0)
    let b = bag(TsallisEntropy(2) => 0.4, TsallisEntropy(3) => 0.6)
        @test length(b) == 2
        @test b[AbstractQAtlas._as_key(TsallisEntropy(2))] == 0.4
    end
    # ...and it does not collide with a Renyi entropy of the same order: the TYPE is
    # still half the key.
    @test AbstractQAtlas._as_key(TsallisEntropy(2)) !=
        AbstractQAtlas._as_key(RenyiEntropy(2))
    @test length(bag(TsallisEntropy(2) => 0.4, RenyiEntropy(2) => 0.6)) == 2
end

@testset "RenyiEntropy carries its order, and validates it" begin
    # QAtlas defines its own `RenyiEntropy` with this field, so the two same-named
    # types were not the same type and `relations_constraining` returned 3 for ours
    # and 0 for theirs -- the three Renyi relations could never see an atlas value.
    # Adopting the field here is what lets QAtlas drop to this one.
    @test RenyiEntropy(2).α == 2.0
    @test RenyiEntropy(2) != RenyiEntropy(3)
    @test_throws ArgumentError RenyiEntropy(0)
    @test_throws ArgumentError RenyiEntropy(-1)
    # α = 1 is the von Neumann limit, refused rather than aliased
    @test_throws ArgumentError RenyiEntropy(1)
end

@testset "a one-parameter family keys under its order, not just its type" begin
    # MEASURED before `OrderSupport` existed:
    #   _as_key(RenyiEntropy(2)) == _as_key(RenyiEntropy(3))   ->  true
    #   bag(RenyiEntropy(2) => 0.5, RenyiEntropy(3) => 0.7)    ->  length 1, 0.7 only
    # A plain field does not distinguish bag keys, because `_as_key` builds the key
    # from `typeof(v)` and `typeof` erases a non-parametric struct's field. The order
    # therefore lives in the support slot, which already existed for "same quantity,
    # different instance" (`RegionSupport`).
    k2 = AbstractQAtlas._as_key(RenyiEntropy(2))
    k3 = AbstractQAtlas._as_key(RenyiEntropy(3))
    @test k2 != k3
    @test k2.type === RenyiEntropy && k3.type === RenyiEntropy    # same type...
    @test k2.support == OrderSupport(2.0)                          # ...different support
    @test k2 == AbstractQAtlas._as_key(RenyiEntropy(2))            # and it is value-based
    @test hash(k2) == hash(AbstractQAtlas._as_key(RenyiEntropy(2)))

    b = bag(RenyiEntropy(2) => 0.5, RenyiEntropy(3) => 0.7)
    @test length(b) == 2
    @test b[k2] == 0.5 && b[k3] == 0.7

    # quantities whose type IS their whole identity still key globally
    @test AbstractQAtlas._as_key(VonNeumannEntropy()).support isa Global
    @test variable_support(VonNeumannEntropy()) isa Global
    @test variable_support(RenyiEntropy(2)) isa OrderSupport
end

@testset "bag refuses a duplicate key instead of overwriting" begin
    # The general form of the bug above: a bag that quietly holds fewer values than
    # it was handed. Now an error wherever the collision comes from.
    @test_throws ErrorException bag(RenyiEntropy(2) => 0.5, RenyiEntropy(2) => 0.7)
    @test_throws ErrorException bag(FreeEnergy => 1.0, FreeEnergy => 2.0)
    # ...and a genuine pair of distinct keys is unaffected
    @test length(bag(FreeEnergy => 1.0, Energy => 2.0)) == 2
end

@testset "sector gaps, band velocities and the T=0 entropy are their own slots" begin
    # Vocabulary ported up from QAtlas (QAtlasHub/QAtlas.jl#807): the families
    # already lived here while the members lived downstream, which put the vague
    # name (MassGap) in the base package and the sector-resolved ones in the atlas.
    @test ChargeGap() isa AbstractGap
    @test SpinGap() isa AbstractGap
    @test FermiVelocity() isa AbstractVelocity
    @test LuttingerVelocity() isa AbstractVelocity
    @test ResidualEntropy() isa AbstractThermalPotential
    @test LogarithmicNegativity() isa AbstractEntanglementMeasure
    @test PageEntropy() isa AbstractEntanglementMeasure

    # The alias is one quantity under two names, not a second slot.
    @test SpinWaveVelocity === LuttingerVelocity

    # ...and the sector gaps are NOT MassGap. Asserted on the identity a bag keys
    # on, because that is what a relation sees: three gaps must occupy three slots.
    # A future "just alias ChargeGap to MassGap" would collapse this to 1 and now
    # fails loudly (bag refuses duplicate keys) instead of silently dropping values.
    @test length(bag(MassGap => 0.0, ChargeGap => 2.0, SpinGap => 0.0)) == 3
    @test length(bag(FermiVelocity => 1.0, LuttingerVelocity => 1.5)) == 2
    @test length(bag(ThermalEntropy => 0.4, ResidualEntropy => 0.3230659669)) == 2
    @test variable_support(ChargeGap()) isa Global
end
