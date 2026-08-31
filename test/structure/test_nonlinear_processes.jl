using AbstractQAtlas
using AbstractQAtlas:
    CurrentCorrelation,
    DynamicalConductivity,
    DynamicalCorrelation,
    frequency_arguments,
    index_spaces,
    response_order,
    tensor_rank
using Test

const PROCESSES_1F = (
    HarmonicGeneration(1),
    HarmonicGeneration(2),
    HarmonicGeneration(3),
    OpticalRectification(),
)
const PROCESSES_2F = (SumFrequencyGeneration(), DifferenceFrequencyGeneration())

@testset "nonlinear processes: a cut supplies exactly as many frequencies as the order" begin
    # The structural law of this layer. A cut that supplied the wrong number of arguments
    # would be silently accepted by any consumer that splats it into χ⁽ⁿ⁾.
    for p in PROCESSES_1F
        @test length(process_frequencies(p, 0.37)) == response_order(p)
    end
    for p in PROCESSES_2F
        @test length(process_frequencies(p, 0.37, 1.4)) == response_order(p)
    end
end

@testset "nonlinear processes: the cut agrees with the quantity on how many frequencies" begin
    # Both sides of the edge must match: the order-n DynamicalSusceptibility carries n
    # frequency arguments, and an order-n process must hand it exactly that many. A cut
    # defined against a frequency-collapsed quantity is the defect this pins against.
    @test frequency_arguments(DynamicalSusceptibility(:x, :y)) ==
        response_order(HarmonicGeneration(1))
    @test frequency_arguments(DynamicalSusceptibility(:x, :y, :z)) ==
        response_order(HarmonicGeneration(2)) ==
        response_order(OpticalRectification()) ==
        response_order(SumFrequencyGeneration()) ==
        response_order(DifferenceFrequencyGeneration())
    @test frequency_arguments(DynamicalSusceptibility(:x, :y, :z, :x)) ==
        response_order(HarmonicGeneration(3))
end

@testset "nonlinear processes: the named cuts are special cases of one another" begin
    # Independent expectations: each identity is derivable from the definitions, so an
    # implementation that got one cut wrong cannot satisfy all of them.
    w, w2 = 0.37, 1.4
    @test process_frequencies(HarmonicGeneration(2), w) ==
        process_frequencies(SumFrequencyGeneration(), w, w)
    @test process_frequencies(OpticalRectification(), w) ==
        process_frequencies(DifferenceFrequencyGeneration(), w, w)
    @test process_frequencies(OpticalRectification(), w) ==
        process_frequencies(SumFrequencyGeneration(), w, -w)
    @test process_frequencies(HarmonicGeneration(1), w) == (w,)
    @test process_frequencies(SumFrequencyGeneration(), w, w2) == (w, w2)
    @test process_frequencies(DifferenceFrequencyGeneration(), w, w2) == (w, -w2)
end

@testset "nonlinear processes: emitted frequency is the sum, not the drive" begin
    w, w2 = 0.37, 1.4
    @test emitted_frequency(HarmonicGeneration(2), w) ≈ 2w
    @test emitted_frequency(HarmonicGeneration(3), w) ≈ 3w
    @test emitted_frequency(SumFrequencyGeneration(), w, w2) ≈ w + w2
    @test emitted_frequency(DifferenceFrequencyGeneration(), w, w2) ≈ w - w2
    # the case the two-function split exists for: nonzero arguments, zero emission
    @test emitted_frequency(OpticalRectification(), w) ≈ 0
    @test all(!iszero, process_frequencies(OpticalRectification(), w))
end

@testset "nonlinear processes: the cut respects intrinsic permutation symmetry" begin
    # χ⁽²⁾(ω₁,ω₂) = χ⁽²⁾(ω₂,ω₁), so a second-order cut cannot depend on the order in which
    # its two drives are named — the frequencies must merely permute, and the emission must
    # be invariant. Ties this layer to structure/tensor_symmetry.jl rather than restating it.
    @test intrinsic_permutation_symmetric(DynamicalSusceptibility(:x, :y, :z))
    a, b = 0.37, 1.4
    fwd = process_frequencies(SumFrequencyGeneration(), a, b)
    rev = process_frequencies(SumFrequencyGeneration(), b, a)
    @test sort(collect(fwd)) == sort(collect(rev))
    @test emitted_frequency(SumFrequencyGeneration(), a, b) ≈
        emitted_frequency(SumFrequencyGeneration(), b, a)
end

@testset "nonlinear processes: parity selection rule" begin
    # An inversion-odd observable (the electric current) driven by an inversion-odd field:
    # linear response survives, second order is forbidden, third survives. This is why a
    # centrosymmetric crystal has no second-harmonic generation.
    @test !parity_forbidden(1, -1)
    @test parity_forbidden(2, -1)
    @test !parity_forbidden(3, -1)
    # an inversion-even observable is the mirror image
    @test parity_forbidden(1, +1)
    @test !parity_forbidden(2, +1)
    @test parity_forbidden(3, +1)
    # and the zeroth order (the equilibrium value) is never forbidden for an even observable
    @test !parity_forbidden(0, +1)

    @test parity_forbidden(HarmonicGeneration(2), -1)
    @test parity_forbidden(OpticalRectification(), -1)
    @test !parity_forbidden(HarmonicGeneration(3), -1)
    # reading the order off a quantity gives the same answer as passing it directly
    chi2 = DynamicalSusceptibility(:x, :y, :z)
    @test response_order(chi2) == 2
    @test parity_forbidden(chi2, -1) == parity_forbidden(2, -1)
end

@testset "nonlinear processes: guards" begin
    @test_throws ArgumentError HarmonicGeneration(0)
    @test_throws ArgumentError HarmonicGeneration(-2)
    @test_throws ArgumentError parity_forbidden(2, 0)      # parity must be ±1
    @test_throws ArgumentError parity_forbidden(2, 2)
    @test_throws ArgumentError parity_forbidden(-1, 1)     # order must be non-negative
end

# ── generality: the layer is one object, not a list of special cases ──────────────────

# an independent count of the distinct arrangements, by enumeration rather than by the
# multinomial formula degeneracy_factor uses.
function _all_perms(v)
    length(v) <= 1 && return [collect(v)]
    out = Vector{Vector{eltype(v)}}()
    for i in eachindex(v)
        rest = vcat(v[1:(i - 1)], v[(i + 1):end])
        for p in _all_perms(rest)
            push!(out, vcat(v[i], p))
        end
    end
    return out
end
_distinct_arrangements(args) = length(unique(_all_perms(collect(args))))

@testset "wave mixing: every named process is one general object" begin
    # No per-process machinery: each name is a constructor for a multiplicity pattern, so
    # arbitrary order and arbitrary colour count are covered by construction.
    @test HarmonicGeneration(2) === WaveMixing((2,), (0,))
    @test OpticalRectification() === WaveMixing((1,), (1,))
    @test SumFrequencyGeneration() === WaveMixing((1, 1), (0, 0))
    @test DifferenceFrequencyGeneration() === WaveMixing((1, 0), (0, 1))
    @test FourWaveMixing() === WaveMixing((2, 0), (0, 1))
    @test KerrEffect() === WaveMixing((2,), (1,))
    @test CrossPhaseModulation() === WaveMixing((1, 1), (1, 0))

    # high-harmonic generation needs no new type
    @test response_order(HarmonicGeneration(17)) == 17
    @test emitted_frequency(HarmonicGeneration(17), 0.1) ≈ 1.7
    @test n_drives(HarmonicGeneration(17)) == 1
    # neither does an arbitrary five-wave, three-colour process
    exotic = WaveMixing((2, 1, 0), (0, 1, 1))
    @test response_order(exotic) == 5
    @test n_drives(exotic) == 3
    @test process_frequencies(exotic, 1, 3, 7) == (1, 1, 3, -3, -7)
    @test emitted_frequency(exotic, 1, 3, 7) == -5
end

@testset "wave mixing: the third-order named processes" begin
    a, b = 0.4, 1.1
    @test process_frequencies(KerrEffect(), a) == (a, a, -a)
    @test emitted_frequency(KerrEffect(), a) ≈ a               # emits back at the drive
    @test process_frequencies(FourWaveMixing(), a, b) == (a, a, -b)
    @test emitted_frequency(FourWaveMixing(), a, b) ≈ 2a - b
    @test process_frequencies(CrossPhaseModulation(), a, b) == (a, -a, b)
    @test emitted_frequency(CrossPhaseModulation(), a, b) ≈ b  # pump contributes no net ω
    @test all(
        p -> response_order(p) == 3,
        (KerrEffect(), FourWaveMixing(), CrossPhaseModulation()),
    )
end

@testset "wave mixing: degeneracy factor counts distinct arrangements" begin
    # Checked against an independent enumeration, on GENERIC (non-colliding) drives.
    for (p, ws) in (
        (HarmonicGeneration(2), (1,)),
        (HarmonicGeneration(3), (1,)),
        (OpticalRectification(), (1,)),
        (SumFrequencyGeneration(), (1, 3)),
        (DifferenceFrequencyGeneration(), (1, 3)),
        (KerrEffect(), (1,)),
        (FourWaveMixing(), (1, 3)),
        (CrossPhaseModulation(), (1, 3)),
        (WaveMixing((2, 1, 0), (0, 1, 1)), (1, 3, 7)),
    )
        @test degeneracy_factor(p) == _distinct_arrangements(process_frequencies(p, ws...))
    end
    # the textbook values
    @test degeneracy_factor(HarmonicGeneration(2)) == 1
    @test degeneracy_factor(SumFrequencyGeneration()) == 2
    @test degeneracy_factor(OpticalRectification()) == 2
    @test degeneracy_factor(HarmonicGeneration(3)) == 1
    @test degeneracy_factor(KerrEffect()) == 3
end

@testset "wave mixing: degeneracy belongs to the process, not to where it is evaluated" begin
    # SFG evaluated at ω₂ = ω₁ produces the same ARGUMENTS as second-harmonic generation but
    # is still a different process with a different degeneracy. Conflating the two is exactly
    # the prefactor error this layer exists to prevent, so it is pinned rather than assumed.
    @test process_frequencies(SumFrequencyGeneration(), 1, 1) ==
        process_frequencies(HarmonicGeneration(2), 1)
    @test degeneracy_factor(SumFrequencyGeneration()) !=
        degeneracy_factor(HarmonicGeneration(2))
    @test SumFrequencyGeneration() !== HarmonicGeneration(2)
end

@testset "wave mixing: arity and construction guards" begin
    @test_throws ArgumentError process_frequencies(SumFrequencyGeneration(), 0.5)
    @test_throws ArgumentError process_frequencies(HarmonicGeneration(2), 0.5, 1.0)
    @test_throws ArgumentError WaveMixing((0,), (0,))            # a drive that never enters
    @test_throws ArgumentError WaveMixing((1, 0), (0, 0))        # ... in second position
    @test_throws ArgumentError WaveMixing((-1,), (2,))           # negative multiplicity
end

@testset "response kernel: the time-domain side of the susceptibility" begin
    # The object a real-time method computes had no name; declaring it closes a one-sided
    # Fourier edge (the spin and current CORRELATION channels already had both sides).
    chi = DynamicalSusceptibility(:x, :y, :z)
    ker = ResponseKernel(:x, :y, :z)
    @test AbstractQAtlas.fourier_conjugate_quantity(chi) === typeof(ker)
    @test AbstractQAtlas.fourier_conjugate_quantity(ker) === typeof(chi)
    @test AbstractQAtlas.fourier_pair(chi, ker)
    @test AbstractQAtlas.fourier_pair(ker, chi)                  # conjugacy is symmetric
    @test representation(ker) == (RealSpace(), TimeDomain())
    @test representation(chi) == (MomentumSpace(), FrequencyDomain())
    # order is preserved across the edge — both sides carry n arguments for order n
    @test response_order(ker) == response_order(chi) == 2
    @test frequency_arguments(ker) == frequency_arguments(chi) == 2

    sig = DynamicalConductivity(:x, :y, :z)
    jker = CurrentResponseKernel(:x, :y, :z)
    @test AbstractQAtlas.fourier_pair(sig, jker)
    @test AbstractQAtlas.fourier_pair(jker, sig)
end

@testset "response kernel: causal ordering breaks the symmetry its transform has" begin
    # χ̄⁽ⁿ⁾(ω⃗) is symmetrised; the kernel is supported on 0 ≤ t̄₁ ≤ ⋯ ≤ t̄ₙ and therefore is
    # not. Declaring the kernel permutation-symmetric would be the one-sided error.
    @test intrinsic_permutation_symmetric(DynamicalSusceptibility(:x, :y, :z))
    @test !intrinsic_permutation_symmetric(ResponseKernel(:x, :y, :z))

    @test causally_ordered(0.0, 1.0, 2.0)
    @test causally_ordered(0.5, 0.5)             # coincident times are inside the closure
    @test !causally_ordered(2.0, 1.0)            # out of order
    @test !causally_ordered(-0.1, 1.0)           # before the field acted
    @test causally_ordered(3.0)
    @test_throws ArgumentError causally_ordered()
end

# ── regressions found in review ────────────────────────────────────────────────────────

@testset "degeneracy_factor is exact, and fails loudly rather than wrapping" begin
    # It used to accumulate the multinomial in Int64. Base checks each individual binomial,
    # but their PRODUCT wraps from order ~36 up — silently, and for some shapes NEGATIVE.
    # A negative count is impossible, and nothing downstream would have caught it.
    exact(p) = begin
        n = response_order(p)
        d = big(1)
        for m in Iterators.flatten((p.plus, p.minus))
            d *= binomial(big(n), big(m))
            n -= m
        end
        d
    end
    for p in (
        HarmonicGeneration(17),
        WaveMixing((2, 1, 0), (0, 1, 1)),
        WaveMixing((8, 8), (8, 0)),            # order 24 — well past the old wrap point's shape
        WaveMixing((10, 10), (5, 5)),          # order 30, still inside Int
    )
        @test degeneracy_factor(p) == exact(p)
        @test degeneracy_factor(p) > 0
    end
    # beyond Int the answer is refused, not wrapped, and the true value is in the message
    big_one = WaveMixing((6, 9), (10, 11))     # order 36; used to return −8617483035347986816
    @test exact(big_one) > typemax(Int)
    @test_throws OverflowError degeneracy_factor(big_one)
    @test_throws "9829261038361564800" degeneracy_factor(big_one)
end

@testset "drive frequencies of mixed numeric type are promoted, not rejected" begin
    # The WaveMixing docstring says the Pockels effect is this process at ω₂ = 0. Written
    # that way — an Int literal beside a Float — it used to die inside Base's tuple
    # construction with a raw TypeError, because NTuple{n} demands one element type.
    @test process_frequencies(SumFrequencyGeneration(), 0.5, 0) === (0.5, 0.0)
    @test emitted_frequency(SumFrequencyGeneration(), 0.5, 0) ≈ 0.5      # Pockels
    @test process_frequencies(SumFrequencyGeneration(), 0.5, 1//2) === (0.5, 0.5)
    @test process_frequencies(FourWaveMixing(), 1, 2.5) === (1.0, 1.0, -2.5)
    # homogeneous input keeps its type
    @test process_frequencies(SumFrequencyGeneration(), 1, 2) === (1, 2)
end

@testset "causally_ordered refuses non-finite times instead of answering" begin
    # NaN used to return `false` — indistinguishable from "genuinely out of order" — and Inf
    # `true`, i.e. an infinite delay reported as inside the support.
    @test_throws ArgumentError causally_ordered(NaN, 1.0)
    @test_throws ArgumentError causally_ordered(0.0, NaN, 2.0)
    @test_throws ArgumentError causally_ordered(0.0, Inf)
    @test_throws "must be finite" causally_ordered(0.0, 1.0, Inf)
    # times are real; a complex argument is a MethodError on the signature, not a comparison
    @test_throws MethodError causally_ordered(1.0 + 0im, 2.0 + 0im)
end

@testset "the Fourier edges are visible to the reflection sweep, not just to this file" begin
    # test/core/test_invariants.jl reaches leaves as bare UnionAlls and skips any type whose
    # fourier_conjugate_quantity returns nothing. Declaring only `::Type{X{I}} where {I}` made
    # all four of these invisible to it — the test that exists to catch a one-sided
    # declaration ran zero assertions against them while reporting green.
    fcq = AbstractQAtlas.fourier_conjugate_quantity
    for (fam, conj) in (
        (DynamicalSusceptibility, ResponseKernel),
        (ResponseKernel, DynamicalSusceptibility),
        (DynamicalConductivity, CurrentResponseKernel),
        (CurrentResponseKernel, DynamicalConductivity),
    )
        @test fcq(fam) === conj                       # the bare form the sweep actually calls
        @test fcq(fam) !== nothing                    # ⇒ the sweep does not `continue` past it
    end
end

@testset "the kernels carry the Kubo edge their docstrings claim" begin
    # ResponseKernel's docstring says it IS the n-fold nested commutator; without an edge the
    # graph showed it reachable only by detouring through the frequency-domain object.
    so = AbstractQAtlas.spectral_origin
    @test so(ResponseKernel(:x, :y, :z)).from === DynamicalCorrelation{(:x, :y, :z)}
    @test so(ResponseKernel(:x, :y, :z)).via === :kubo
    @test so(CurrentResponseKernel(:x, :y, :z)).from === CurrentCorrelation{(:x, :y, :z)}
    @test so(ResponseKernel).from === DynamicalCorrelation          # index-erased form too
    @test so(CurrentResponseKernel).from === CurrentCorrelation
    # order-faithful: the edge preserves the response order, as the χ edge does
    @test response_order(so(ResponseKernel(:x, :y, :z)).from) ==
        response_order(ResponseKernel(:x, :y, :z))
end

@testset "the two kernels live in different index spaces" begin
    # They are near-identical copies of one another; swapping SpinAxis for SpatialDirection is
    # a copy-paste away and length-only trait checks cannot see it.
    @test index_spaces(ResponseKernel(:x, :y, :z)) == (SpinAxis(), SpinAxis(), SpinAxis())
    @test index_spaces(CurrentResponseKernel(:x, :y, :z)) ==
        (SpatialDirection(), SpatialDirection(), SpatialDirection())
    @test index_spaces(ResponseKernel(:x, :y, :z)) !=
        index_spaces(CurrentResponseKernel(:x, :y, :z))
    # and the current kernel's traits are asserted directly, not only via its transform
    @test response_order(CurrentResponseKernel(:x, :y, :z)) == 2
    @test frequency_arguments(CurrentResponseKernel(:x, :y, :z)) == 2
    @test tensor_rank(CurrentResponseKernel(:x, :y, :z)) == 3
    @test_throws ErrorException ResponseKernel(:x)            # ≥2 indices
    @test_throws ErrorException CurrentResponseKernel(:x)
end

@testset "HarmonicGeneration validates its own argument, not by accident" begin
    # WaveMixing's own guards also reject (0,)/(−2,), for unrelated reasons, so asserting
    # only "an ArgumentError happened" would still pass if HarmonicGeneration's check were
    # deleted. Pin the message that is actually its own.
    @test_throws "harmonic order must be at least 1" HarmonicGeneration(0)
    @test_throws "harmonic order must be at least 1" HarmonicGeneration(-2)
    @test_throws "never enters" WaveMixing((0,), (0,))
    @test_throws "non-negative" WaveMixing((-2,), (0,))
end
