# structure/criticality.jl — the quantity ⇄ critical-exponent correspondence.
#
# Near a continuous transition each observable has a leading singularity
# governed by ONE of the critical exponents.  That correspondence —
# which exponent belongs to which quantity, and with what power — is a
# generic, model-independent structural fact.  It is declared ONCE per
# quantity here; the actual singular forms and the finite-size-scaling
# combinations are then DERIVED, never restated.
#
# This is the layer the previous `universality_forms.jl` was missing:
# there the knowledge "susceptibility scales as L^{γ/ν}" lived in a
# docstring and γ/ν was a hand-passed number.  Here `Susceptibility ↦ γ`
# is code, and `L^{γ/ν}` falls out of it.

"""
    CriticalScaling(exponent, power)

The reduced-temperature critical law of a quantity: `Q ∼ |t|^{power·e}`
where `t = (T − T_c)/T_c` and `e` is the critical exponent named
`exponent` (a field of a [`CriticalExponents`](@ref) NamedTuple).
`power = +1` for a quantity that *vanishes* at criticality (e.g. the
order parameter, `M ∼ |t|^{+β}`), `power = −1` for one that *diverges*
(e.g. `χ ∼ |t|^{−γ}`, `ξ ∼ |t|^{−ν}`, `C ∼ |t|^{−α}`).
"""
struct CriticalScaling
    exponent::Symbol
    power::Int
end

"""
    critical_scaling(quantity) -> Union{CriticalScaling,Nothing}
    critical_scaling(::Type{<:AbstractQuantity}) -> Union{CriticalScaling,Nothing}

The reduced-temperature critical correspondence of `quantity`: which
exponent governs its `|t|`-singularity, and with what sign — or
`nothing` if the quantity has no reduced-temperature critical law
(e.g. the partition function, or the field-driven δ / distance-driven η
laws handled by [`critical_isotherm`](@ref) / [`correlation_decay`](@ref)).

```julia
critical_scaling(Susceptibility)          # CriticalScaling(:γ, -1) ⇒ χ ∼ |t|^{-γ}
critical_scaling(SpontaneousMagnetization) # CriticalScaling(:β, +1) ⇒ M ∼ |t|^{+β}
```
"""
critical_scaling(q::AbstractQuantity) = critical_scaling(typeof(q))
critical_scaling(::Type{<:AbstractQuantity}) = nothing

critical_scaling(::Type{SpontaneousMagnetization}) = CriticalScaling(:β, +1)
critical_scaling(::Type{<:AbstractSusceptibility}) = CriticalScaling(:γ, -1)
critical_scaling(::Type{SpecificHeat}) = CriticalScaling(:α, -1)
critical_scaling(::Type{CorrelationLength}) = CriticalScaling(:ν, -1)
export critical_scaling, CriticalScaling

"""
    critical_isotherm(::Type{SpontaneousMagnetization}) -> Symbol

The critical-isotherm exponent: exactly at `T_c`, the order parameter
responds to its conjugate field as `M ∼ h^{1/δ}`.  Returns the exponent
symbol `:δ`.  A distinct functional form from the reduced-temperature
laws, hence its own accessor.
"""
critical_isotherm(::Type{SpontaneousMagnetization}) = :δ
critical_isotherm(q::AbstractQuantity) = critical_isotherm(typeof(q))
export critical_isotherm

"""
    correlation_decay(::Type{<:AbstractTwoPointCorrelation}) -> Symbol

The anomalous-dimension exponent: at `T_c` the two-point function decays
as `G(r) ∼ r^{−(d−2+η)}`.  Returns the exponent symbol `:η`
(the full decay power needs the spatial dimension `d`, supplied at use).
"""
correlation_decay(::Type{<:AbstractTwoPointCorrelation}) = :η
correlation_decay(q::AbstractQuantity) = correlation_decay(typeof(q))
export correlation_decay

# ─── Derived singular forms (correspondence-driven, nothing hand-passed) ─

"""
    singular_form(quantity, t; exponents::NamedTuple) -> Real

The leading singular form of `quantity` at reduced temperature `t`,
`Q ∼ |t|^{power·e}`, with the exponent looked up from the
[`critical_scaling`](@ref) correspondence and its value taken from
`exponents`.  The correspondence — not the caller — decides which
exponent and which sign:

```julia
exps = (α=0//1, β=1//8, γ=7//4, δ=15//1, ν=1//1, η=1//4)
singular_form(SpontaneousMagnetization(), -0.01; exponents=exps)  # |t|^{+1/8}
singular_form(SusceptibilityZZ(), 0.01; exponents=exps)           # |t|^{-7/4}
```

Throws for a quantity with no reduced-temperature critical law.
"""
function singular_form(q::AbstractQuantity, t; exponents::NamedTuple)
    cs = critical_scaling(q)
    cs === nothing && error(
        "singular_form: $(typeof(q)) has no reduced-temperature critical law " *
        "(critical_scaling returned nothing)",
    )
    return abs(t)^(cs.power * exponents[cs.exponent])
end
export singular_form

"""
    fss_size_exponent(quantity; exponents::NamedTuple) -> Real

The exponent of the linear size `L` in the finite-size scaling of
`quantity` at criticality: `Q(T_c, L) ∼ L^{fss_size_exponent}`.  Derived
from the correspondence — a divergent quantity `Q ∼ |t|^{−x}` scales as
`L^{+x/ν}`, a vanishing one `Q ∼ |t|^{+x}` as `L^{−x/ν}` — so the size
exponent is always `−(power·e)/ν`:

```julia
fss_size_exponent(SusceptibilityZZ(); exponents=exps)          # +γ/ν  (= 7/4)
fss_size_exponent(SpontaneousMagnetization(); exponents=exps)  # −β/ν  (= −1/8)
fss_size_exponent(CorrelationLength(); exponents=exps)         # +1    (ξ ∼ L)
```
"""
function fss_size_exponent(q::AbstractQuantity; exponents::NamedTuple)
    cs = critical_scaling(q)
    cs === nothing &&
        error("fss_size_exponent: $(typeof(q)) has no reduced-temperature critical law")
    return -(cs.power * exponents[cs.exponent]) / exponents.ν
end
export fss_size_exponent

"""
    fss_peak(quantity, L; exponents::NamedTuple) -> Real

Finite-size scaling of `quantity`'s critical peak with linear size `L`,
`Q(T_c, L) ∝ L^{fss_size_exponent}` — the exponent derived from the
[`critical_scaling`](@ref) correspondence, not passed by hand.

```julia
fss_peak(SusceptibilityZZ(), 64; exponents=exps)   # ∝ 64^{γ/ν}
```
"""
function fss_peak(q::AbstractQuantity, L; exponents::NamedTuple)
    return L^fss_size_exponent(q; exponents=exponents)
end
export fss_peak

"""
    collapse_coordinates(quantity, T, L, Tc; exponents::NamedTuple) -> (x, scale)

The finite-size-scaling data-collapse transform for `quantity`, with
both exponent combinations derived from the correspondence.  For an
observable obeying `O(T, L) = L^{−ρ} f((T − T_c)·L^{1/ν})` (with
`ρ = −fss_size_exponent`), plotting `O·scale` against `x` collapses all
sizes onto the universal curve `f`:

- `x     = (T − Tc)·L^{1/ν}`,
- `scale = L^{−fss_size_exponent}` = `L^{ρ}`.

At `T = Tc`, `x = 0` for every `L` (the collapse pivot).  The residual
spread of the collapsed data across sizes is the quantitative
universality test — and the exponents used are exactly the atlas's,
never hand-typed.
"""
function collapse_coordinates(q::AbstractQuantity, T, L, Tc; exponents::NamedTuple)
    ρ = -fss_size_exponent(q; exponents=exponents)
    return (x=(T - Tc) * L^(1 / exponents.ν), scale=L^ρ)
end
export collapse_coordinates
