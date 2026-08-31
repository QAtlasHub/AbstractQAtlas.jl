# structure/nonlinear_processes.jl — which frequency arguments a nonlinear process IS.
#
# `χ⁽ⁿ⁾(ω₁, …, ωₙ)` is one object; "second-harmonic generation", "optical rectification",
# "four-wave mixing" are not separate quantities but *cuts* through it — particular choices
# of its frequency arguments.  Every consumer that plots or extracts one of them has to
# encode that choice, and a consumer that encodes it locally ("the second harmonic is the
# point `ω = 2ω₁`") has silently forked the convention.  This layer states each cut once, so
# a downstream extraction asks rather than re-derives.
#
# The general object is [`WaveMixing`](@ref): `M` drive frequencies, each entering the
# response a given number of times with `+ω` and a given number of times with `−ω`.  Every
# named process is a constructor for one — there is no separate machinery per process, so
# "is high-harmonic generation supported?" is not a question about coverage.
#
# Three functions carry it.  `process_frequencies` gives the arguments handed to `χ⁽ⁿ⁾`;
# `emitted_frequency` gives their sum — the frequency at which the induced response actually
# appears, which is what a measurement indexes by; `degeneracy_factor` counts the distinct
# arrangements of those arguments, which is the bookkeeping between the KERNEL and what a
# monochromatic drive produces.  Keeping all three explicit is the point: they differ
# (rectification emits at 0 from nonzero arguments and has degeneracy 2), and conflating any
# two of them is the mistake this layer exists to prevent.
#
# The cuts compose with, and are independent of, the intrinsic permutation symmetry in
# `structure/tensor_symmetry.jl`: that permutes (field-index, frequency) PAIRS, this selects
# the frequencies.

"""
    AbstractNonlinearProcess

A named nonlinear process, i.e. a choice of frequency arguments for `χ⁽ⁿ⁾(ω₁, …, ωₙ)`.
[`WaveMixing`](@ref) is the general (and currently only) concrete form.

A process is not a quantity — it is a cut through one.  The quantity is the order-`n`
[`DynamicalSusceptibility`](@ref) (or [`DynamicalConductivity`](@ref) in the current
channel); the process says which of its frequency arguments to evaluate.
"""
abstract type AbstractNonlinearProcess end
export AbstractNonlinearProcess

"""
    WaveMixing(plus, minus)

A general wave-mixing process: `M = length(plus)` drive frequencies, where drive `k` enters
the response `plus[k]` times as `+ω_k` and `minus[k]` times as `−ω_k`.  The total response
order is `sum(plus) + sum(minus)` and the emitted frequency is `Σ (plus[k] − minus[k]) ω_k`.

Every named process below is a constructor for one of these, so the layer covers arbitrary
order and arbitrary numbers of colours without new machinery:

| process | `WaveMixing` | arguments | emits |
|---|---|---|---|
| [`HarmonicGeneration`](@ref)`(q)` | `((q,), (0,))` | `(ω, …, ω)` | `qω` |
| [`OpticalRectification`](@ref) | `((1,), (1,))` | `(ω, −ω)` | `0` |
| [`SumFrequencyGeneration`](@ref) | `((1,1), (0,0))` | `(ω₁, ω₂)` | `ω₁+ω₂` |
| [`DifferenceFrequencyGeneration`](@ref) | `((1,0), (0,1))` | `(ω₁, −ω₂)` | `ω₁−ω₂` |
| [`FourWaveMixing`](@ref) | `((2,0), (0,1))` | `(ω₁, ω₁, −ω₂)` | `2ω₁−ω₂` |
| [`KerrEffect`](@ref) | `((2,), (1,))` | `(ω, ω, −ω)` | `ω` |
| [`CrossPhaseModulation`](@ref) | `((1,1), (1,0))` | `(ω₁, −ω₁, ω₂)` | `ω₂` |

Processes that are a *special evaluation point* rather than a distinct multiplicity pattern
do not need their own name: the Pockels effect is
`SumFrequencyGeneration()` at `ω₂ = 0`, and electric-field-induced second-harmonic
generation is `FourWaveMixing()`-shaped with one drive held at zero frequency.

A drive that never enters is rejected — it would make the arity of
[`process_frequencies`](@ref) disagree with the process being described.
"""
struct WaveMixing{M} <: AbstractNonlinearProcess
    plus::NTuple{M,Int}
    minus::NTuple{M,Int}
    function WaveMixing(plus::NTuple{M,Int}, minus::NTuple{M,Int}) where {M}
        M >= 1 || throw(ArgumentError("a wave-mixing process needs at least one drive"))
        (all(>=(0), plus) && all(>=(0), minus)) || throw(
            ArgumentError("multiplicities must be non-negative, got +$plus and −$minus")
        )
        idle = findfirst(k -> plus[k] + minus[k] == 0, 1:M)
        idle === nothing ||
            throw(ArgumentError("drive $idle never enters the process; drop it"))
        return new{M}(plus, minus)
    end
end
function WaveMixing(plus::NTuple{M,Integer}, minus::NTuple{M,Integer}) where {M}
    return WaveMixing(Int.(plus), Int.(minus))
end
export WaveMixing

"""
    n_drives(process) -> Int

How many independent drive frequencies the process consumes — the number of arguments
[`process_frequencies`](@ref) and [`emitted_frequency`](@ref) expect.
"""
n_drives(::WaveMixing{M}) where {M} = M
export n_drives

response_order(p::WaveMixing) = sum(p.plus) + sum(p.minus)

"""
    HarmonicGeneration(q)

`q`-th harmonic generation: all `q` arguments equal, `χ⁽q⁾(ω, …, ω)`, emitting at `qω`.
`q = 2` is second-harmonic generation ([Franken1961](@cite), the first observation of a
nonlinear optical process); large `q` is high-harmonic generation, which needs no separate
type.
"""
function HarmonicGeneration(order::Integer)
    order >= 1 || throw(ArgumentError("harmonic order must be at least 1, got $order"))
    return WaveMixing((Int(order),), (0,))
end
export HarmonicGeneration

"""
    OpticalRectification()

Optical rectification: `χ⁽²⁾(ω, −ω)`, emitting at **zero** frequency — a static response
induced by an oscillating field ([Bass1962](@cite)).  The photovoltaic / photogalvanic
effect is its current-channel counterpart.

The distinction this layer exists for is visible here: the frequency arguments are `(ω, −ω)`
and nonzero, the emitted frequency is `0`, and the degeneracy factor is `2`.  All three are
different numbers.
"""
OpticalRectification() = WaveMixing((1,), (1,))
export OpticalRectification

"""
    SumFrequencyGeneration()

Sum-frequency generation: `χ⁽²⁾(ω₁, ω₂)`, emitting at `ω₁ + ω₂`.  The generic second-order
two-colour cut — [`HarmonicGeneration`](@ref)`(2)` is its `ω₂ = ω₁` case, and the Pockels
effect is its `ω₂ = 0` case.
"""
SumFrequencyGeneration() = WaveMixing((1, 1), (0, 0))
export SumFrequencyGeneration

"""
    DifferenceFrequencyGeneration()

Difference-frequency generation: `χ⁽²⁾(ω₁, −ω₂)`, emitting at `ω₁ − ω₂`.
[`OpticalRectification`](@ref) is its `ω₂ = ω₁` case.
"""
DifferenceFrequencyGeneration() = WaveMixing((1, 0), (0, 1))
export DifferenceFrequencyGeneration

"""
    FourWaveMixing()

Four-wave mixing: `χ⁽³⁾(ω₁, ω₁, −ω₂)`, emitting at `2ω₁ − ω₂` — the third-order
two-colour process behind coherent anti-Stokes Raman scattering and phase conjugation
([ArmstrongBloembergen1962](@cite) for the permutation bookkeeping it obeys).
"""
FourWaveMixing() = WaveMixing((2, 0), (0, 1))
export FourWaveMixing

"""
    KerrEffect()

The optical Kerr effect / self-phase modulation: `χ⁽³⁾(ω, ω, −ω)`, emitting back at `ω` —
an intensity-dependent refractive index.  Single-colour, third order, degeneracy factor `3`.
"""
KerrEffect() = WaveMixing((2,), (1,))
export KerrEffect

"""
    CrossPhaseModulation()

Cross-phase modulation: `χ⁽³⁾(ω₁, −ω₁, ω₂)`, emitting at `ω₂` — a pump at `ω₁` modulating a
probe at `ω₂`.  The pump enters twice with opposite signs, so it contributes no net
frequency; this is the pump–probe shape.
"""
CrossPhaseModulation() = WaveMixing((1, 1), (1, 0))
export CrossPhaseModulation

"""
    process_frequencies(process, ω...) -> NTuple

The frequency arguments `(ω₁, …, ωₙ)` at which `χ⁽ⁿ⁾` is to be evaluated for `process`,
given the `n_drives(process)` driving frequencies actually applied.  The result always has
`response_order(process)` entries.

```julia
process_frequencies(HarmonicGeneration(3), 0.5)          # (0.5, 0.5, 0.5)
process_frequencies(OpticalRectification(), 0.5)         # (0.5, -0.5)
process_frequencies(SumFrequencyGeneration(), 0.5, 1.2)  # (0.5, 1.2)
process_frequencies(FourWaveMixing(), 0.5, 1.2)          # (0.5, 0.5, -1.2)
```
"""
function process_frequencies(p::WaveMixing{M}, ws::Vararg{Any,N}) where {M,N}
    N == M || throw(
        ArgumentError("this process has $M drive(s) but got $N frequency argument(s)")
    )
    # `NTuple{n}` demands ONE element type, so mixing an Int with a Float — which the
    # `WaveMixing` docstring invites, the Pockels effect being this process at `ω₂ = 0` —
    # otherwise died inside Base's tuple construction with a raw TypeError. Promote first:
    # mixed numeric types are ordinary in every other Julia numeric API.
    fs = promote(ws...)
    args = Iterators.flatten(
        Iterators.flatten((
            Iterators.repeated(fs[k], p.plus[k]), Iterators.repeated(-fs[k], p.minus[k])
        )) for k in 1:M
    )
    return NTuple{response_order(p)}(args)
end
export process_frequencies

"""
    emitted_frequency(process, ω...) -> Number

The frequency at which the induced response appears: the SUM of
[`process_frequencies`](@ref), because the order-`n` response at total frequency `ω`
collects the arguments obeying `ω₁ + ⋯ + ωₙ = ω`.

This is the frequency a measurement indexes by, and it is NOT the driving frequency: optical
rectification emits at `0` while being driven at `ω`, and cross-phase modulation emits at the
probe frequency while being driven by a pump as well.

```julia
emitted_frequency(HarmonicGeneration(2), 0.5)      # 1.0
emitted_frequency(OpticalRectification(), 0.5)     # 0.0
emitted_frequency(FourWaveMixing(), 1, 3)          # -1    ( = 2·1 − 3 )
```
"""
emitted_frequency(p::AbstractNonlinearProcess, w...) = sum(process_frequencies(p, w...))
export emitted_frequency

"""
    degeneracy_factor(process) -> Int

The number of **distinct arrangements** of the process's frequency arguments — the
multinomial `n! / ∏ mᵢ!` over the multiplicities.

This is the bookkeeping between the response *kernel* and what a monochromatic drive
produces: the frequency integral defining the order-`n` response sums over every
arrangement of the applied frequencies that meets the emission condition, so a process whose
arguments are all distinct collects more terms than a fully degenerate one.  Second-harmonic
generation has degeneracy `1` while sum-frequency generation has `2`, which is why the two
carry different prefactors in every textbook table and why the difference is worth stating
once rather than per consumer.

```julia
degeneracy_factor(HarmonicGeneration(2))       # 1  — (ω, ω)
degeneracy_factor(SumFrequencyGeneration())    # 2  — (ω₁, ω₂)
degeneracy_factor(OpticalRectification())      # 2  — (ω, −ω)
degeneracy_factor(KerrEffect())                # 3  — (ω, ω, −ω)
```
"""
function degeneracy_factor(p::WaveMixing)
    # Accumulated in `BigInt`, deliberately. Base's `binomial` checks each individual
    # coefficient, but the running PRODUCT is the full multinomial and wrapped `Int64`
    # silently from order ~36 upward — returning, among other things, NEGATIVE counts. A
    # wrong-but-plausible integer out of a library used as an oracle is worse than an error,
    # so the arithmetic is exact and the single narrowing at the end fails loudly, with the
    # true value in the message.
    d = big(1)
    remaining = response_order(p)
    for m in Iterators.flatten((p.plus, p.minus))
        m == 0 && continue
        d *= binomial(big(remaining), big(m))
        remaining -= m
    end
    typemin(Int) <= d <= typemax(Int) || throw(
        OverflowError(
            "degeneracy_factor: this order-$(response_order(p)) process has $d distinct " *
            "argument arrangements, which does not fit in Int",
        ),
    )
    return Int(d)
end
export degeneracy_factor

"""
    parity_forbidden(order, observable_parity) -> Bool
    parity_forbidden(process_or_quantity, observable_parity) -> Bool

Whether a symmetry forces the order-`n` response to vanish identically.

If the system has a symmetry under which the driving field is **odd** (`f → −f`) and the
observable transforms with parity `s = ±1` (`Q → sQ`), then `χ⁽ⁿ⁾ = s(−1)ⁿ χ⁽ⁿ⁾`, so the
response vanishes unless `s(−1)ⁿ = +1`:

  * an **even** observable (`s = +1`) in such a system has no **odd**-order response;
  * an **odd** observable (`s = −1`) has no **even**-order response.

Inversion in a centrosymmetric crystal is the standard instance: the current is
inversion-odd, so `χ⁽²⁾ ≡ 0` and second-harmonic generation is forbidden.

This is a *selection rule*, not a numeric identity, so it is structure and not an
`AbstractRelation`: it tells a consumer when to EXPECT zero, which is what makes a measured
zero informative rather than vacuous.

```julia
parity_forbidden(2, -1)                          # true  — no χ⁽²⁾ for an odd observable
parity_forbidden(HarmonicGeneration(2), -1)      # true
parity_forbidden(HarmonicGeneration(3), -1)      # false — χ⁽³⁾ survives
```
"""
function parity_forbidden(order::Integer, observable_parity::Integer)
    order >= 0 || throw(ArgumentError("response order must be non-negative, got $order"))
    observable_parity in (-1, 1) ||
        throw(ArgumentError("observable parity must be ±1, got $observable_parity"))
    return observable_parity * (-1)^order != 1
end
function parity_forbidden(p::AbstractNonlinearProcess, observable_parity::Integer)
    return parity_forbidden(response_order(p), observable_parity)
end
function parity_forbidden(q::AbstractQuantity, observable_parity::Integer)
    return parity_forbidden(response_order(q), observable_parity)
end
export parity_forbidden

"""
    causally_ordered(t̄...) -> Bool

Whether the time arguments lie in the support of a retarded response kernel:
`0 ≤ t̄₁ ≤ t̄₂ ≤ ⋯ ≤ t̄ₙ`, with `t̄ᵢ = t − tᵢ` the delay from each field application to the
measurement.

The kernel [`ResponseKernel`](@ref) vanishes outside this region — the fields must act
before the measurement and, in the ordered form, in sequence.  Two consequences a consumer
should not have to rediscover: the kernel is **not** permutation symmetric (the ordering
breaks it), while its Fourier transform is symmetrised into `χ̄⁽ⁿ⁾`; and integrating over
the ordered region is what produces the nested denominators of the frequency-domain
response.
"""
function causally_ordered(ts::Vararg{Real,N}) where {N}
    N >= 1 || throw(ArgumentError("need at least one time argument"))
    # Without this, NaN returns `false` — indistinguishable from a genuinely out-of-order
    # sequence — and Inf returns `true`, telling a caller an infinite delay lies inside the
    # support. Both bury a corrupted time argument under a confident Boolean. The same guard
    # is house style: see the `all(isfinite, ...)` check in relations/interface.jl.
    all(isfinite, ts) ||
        throw(ArgumentError("causally_ordered: time arguments must be finite, got $ts"))
    ts[1] >= 0 || return false
    return all(ts[i] <= ts[i + 1] for i in 1:(N - 1))
end
export causally_ordered
