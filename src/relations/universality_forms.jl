# relations/universality_forms.jl — scaling FORMS.
#
# The standard functional forms of critical scaling, as small documented
# pure functions.  These are maps (what an observable looks like near
# criticality), not identities — hence plain functions rather than
# `AbstractRelation` objects.  They are the vocabulary a finite-size
# scaling analysis consumes: fit data against these forms, extract the
# exponents, then gate the exponents with `relations/scaling.jl`.

"""
    fss_peak_scaling(L; ratio) -> L^ratio

Finite-size scaling of a response peak height with linear system size
`L`, e.g. the susceptibility maximum `χ_max(L) ∝ L^{γ/ν}` (pass
`ratio = γ/ν`) or the specific-heat maximum `c_max(L) ∝ L^{α/ν}`.
"""
fss_peak_scaling(L; ratio) = L^ratio
export fss_peak_scaling

"""
    order_parameter_form(t; β) -> |t|^β

Leading singular form of the order parameter on the ordered side of the
transition, `M(t) ∝ |t|^β` for reduced temperature
`t = (T − T_c)/T_c < 0`.
"""
order_parameter_form(t; β) = abs(t)^β
export order_parameter_form

"""
    correlation_length_form(t; ν) -> |t|^(-ν)

Leading divergence of the correlation length, `ξ(t) ∝ |t|^{−ν}`.
"""
correlation_length_form(t; ν) = abs(t)^(-ν)
export correlation_length_form

"""
    susceptibility_form(t; γ) -> |t|^(-γ)

Leading divergence of the susceptibility, `χ(t) ∝ |t|^{−γ}`.
"""
susceptibility_form(t; γ) = abs(t)^(-γ)
export susceptibility_form

"""
    collapse_coordinates(T, L, Tc; ν, ratio) -> (x, scale)

The standard two-exponent finite-size-scaling collapse transform.  For
an observable obeying `O(T, L) = L^{-ratio} · f((T − T_c) · L^{1/ν})`,
plotting `O · scale` against `x` collapses all system sizes onto the
single universal curve `f`:

- `x     = (T − Tc) · L^(1/ν)` — the scaling variable,
- `scale = L^ratio`            — multiply the observable by this.

Pass `ratio = β/ν` for the order parameter (`M ∝ L^{−β/ν}` at `T_c`) or
`ratio = −γ/ν` for the susceptibility (`χ ∝ L^{+γ/ν}` at `T_c`).  At
`T = Tc` the scaling variable is exactly `x = 0` for every `L` — the
collapse pivot.  The residual spread of the collapsed data across sizes
is the quantitative universality test.
"""
function collapse_coordinates(T, L, Tc; ν, ratio)
    return (x=(T - Tc) * L^(1 / ν), scale=L^ratio)
end
export collapse_coordinates
