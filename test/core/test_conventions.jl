# The convention layer: what a project's numbers are written in, and the
# conversion into what the relations here are written in.
#
# The measurement that justifies the layer is the `region_report` pair below.
# `S = 1.75` on two qubits is legal in bits (the ceiling is 2) and impossible in
# nats (the ceiling is 2 ln 2 = 1.386), so a bits-counting calculation handed to
# this package raw is reported as violating a physical bound it does not violate.

using AbstractQAtlas
using Test: @test, @test_throws, @testset

# A quantity, an axis and a conversion, none of them in `src/`: the claim that
# the layer is a protocol rather than a list of axes is a test, not a docstring.
struct ConventionProbeQuantity <: AbstractQuantity end
struct WholeUnits <: Convention end
struct HalfUnits <: Convention end
AbstractQAtlas.canonical_convention(::Type{ConventionProbeQuantity}) = WholeUnits()
AbstractQAtlas.convert_convention(::WholeUnits, ::HalfUnits, ::Type, v) = 2v

@testset "an axis is declared per quantity, never by supertype" begin
    # The twelve whose ABQ definition contains a logarithm, or is an additive
    # combination of ones that do.
    for Q in (
        VonNeumannEntropy,
        FermionicEntanglementEntropy,
        RenyiEntropy,
        MutualInformation,
        ConditionalEntropy,
        RelativeEntropy,
        MeasurementEntropy,
        MarkovEntropy,
        TripartiteInformation,
        TopologicalEntanglementEntropy,
        LogarithmicNegativity,
        PageEntropy,
    )
        @test canonical_convention(Q) === Nats
    end
    # `(1 - Tr ρ^q)/(q - 1)` has no logarithm, and neither has a concurrence or
    # a tangle built by squaring one. A supertype declaration would give all four
    # a base they do not have.
    for Q in (TsallisEntropy, Concurrence, Tangle, ThreeTangle)
        @test canonical_convention(Q) === nothing
    end
    @test canonical_convention(Temperature) === nothing
end

@testset "a declaration is refused when it claims something it cannot mean" begin
    @test_throws ErrorException conventions(Float64 => Bits)
    @test_throws ErrorException conventions(VonNeumannEntropy => 2)
    @test_throws ErrorException conventions(
        VonNeumannEntropy => Bits, VonNeumannEntropy => Nats
    )
    # Naming a concrete type is a claim about that type, so a type with no axis
    # is an error; naming its supertype is a sweep, and skips it silently.
    @test_throws ErrorException conventions(TsallisEntropy => Bits)
    @test conventions(AbstractEntanglementMeasure => Bits) isa ConventionSet
end

@testset "lookup is most specific first" begin
    cs = conventions(AbstractEntanglementMeasure => Bits, VonNeumannEntropy => Nats)
    @test declared_convention(cs, VonNeumannEntropy) === Nats
    @test declared_convention(cs, RenyiEntropy) === Bits
    @test declared_convention(cs, Temperature) === nothing
end

@testset "conversion" begin
    @test convert_convention(Nats, Bits, VonNeumannEntropy, 1.0) ≈ log(2)
    @test convert_convention(Bits, Nats, VonNeumannEntropy, log(2)) ≈ 1.0
    # Equal conventions never touch the value, so an exact input stays exact.
    v = convert_convention(Nats, Nats, VonNeumannEntropy, 3//2)
    @test v === 3//2
    # Two axes that are not a rescaling of each other have no conversion, and
    # returning the value unchanged there is the failure this refuses.
    @test_throws ErrorException convert_convention(
        Nats, HalfUnits(), VonNeumannEntropy, 1.0
    )
    @test_throws ArgumentError LogBase(1)
    @test_throws ArgumentError LogBase(0)
end

@testset "bag converts on entry, and only what is declared" begin
    cs = conventions(AbstractEntanglementMeasure => Bits)
    b = bag(cs, VonNeumannEntropy() => 3.0, TsallisEntropy(2) => 1.0, Temperature => 2.0)
    @test b[VariableKey(VonNeumannEntropy)] ≈ 3.0 * log(2)
    @test b[VariableKey(TsallisEntropy, OrderSupport(2.0))] == 1.0   # no axis
    @test b[VariableKey(Temperature)] == 2.0                        # not declared
    # An undeclared project is the behaviour this package had before the layer.
    @test bag(VonNeumannEntropy() => 3.0)[VariableKey(VonNeumannEntropy)] == 3.0
    # The plain constructor's refusals still hold through this door.
    @test_throws ErrorException bag(cs, VonNeumannEntropy() => nothing)
    @test_throws ErrorException bag(
        cs, VonNeumannEntropy() => 1.0, VonNeumannEntropy() => 2.0
    )
end

@testset "a new axis needs no change to the package" begin
    cs = conventions(ConventionProbeQuantity => HalfUnits())
    @test bag(cs, ConventionProbeQuantity() => 2.5)[VariableKey(ConventionProbeQuantity)] ==
        5.0
end

@testset "the declaration flips a physical-bound verdict, in both directions" begin
    cs = conventions(AbstractEntanglementMeasure => Bits)
    verdict(b) = all(r -> r.pass, region_report(b; local_dim=2))
    # Two qubits: the ceiling is 2 bits, and 2 ln 2 = 1.386 nats.
    @test !verdict(bag(entanglement_entropy(1, 2) => 1.75))        # false alarm
    @test verdict(bag(cs, entanglement_entropy(1, 2) => 1.75))     # declared away
    # The control: 2.5 bits is over the ceiling in bits too, so declaring the
    # convention must not make it pass.
    @test !verdict(bag(entanglement_entropy(1, 2) => 2.5))
    @test !verdict(bag(cs, entanglement_entropy(1, 2) => 2.5))
end
