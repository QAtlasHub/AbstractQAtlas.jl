# evaluation.jl — the FUNCTIONAL-EVALUATION interface (the scope-line #14 seam).
#
# The supplied-integral relations — Kramers–Kronig (a principal-value Hilbert transform),
# the spectral / structure-factor sum rules (∫A, ∫S) and the f-sum rule (∫ωS) — need a
# TRANSFORM of a response FUNCTION, which this stdlib-only, model-independent package
# deliberately does NOT evaluate: computing a principal-value Hilbert transform or a
# spectral quadrature belongs to the FUNCTIONAL SIBLING (a numerical / ParaLA-based
# package; issues #14, #19).  This file OWNS the generic evaluation verbs — exactly as
# `fetch` owns the value-retrieval verb — with an informative "not implemented" fallback.
# A functional package subtypes `AbstractResponse` for its representation (an (ω, values)
# grid, an analytic pole–residue rep, …) and adds the methods; after that the
# supplied-integral `check`/`residual` become turnkey straight from a measured response —
# `check(SpectralSumRule(); spectral_integral = spectral_moment(A_rep, 0))`, etc.  The
# interface is here; the numerics live there.

"""
    AbstractResponse

Parent type for a representation of a response FUNCTION over frequency (or `(q, ω)`) — an
`(ω, values)` grid, an analytic pole–residue rep, … — that the functional sibling can
transform.  AbstractQAtlas owns only this abstract type and the evaluation verbs
([`principal_value_hilbert`](@ref), [`spectral_moment`](@ref)); the concrete
representations and their methods live in the functional package (the `fetch`-style seam:
the interface here, the numerics there).
"""
abstract type AbstractResponse end
export AbstractResponse

"""
    principal_value_hilbert(response::AbstractResponse, ω) -> Number

The principal-value Hilbert transform `P ∫ f(ω′)/(ω′ − ω) dω′` of `response` at `ω` — the
`pv_real` / `pv_imag` a [`KramersKronigReal`](@ref) / [`KramersKronigImag`](@ref) check
consumes (feed the imaginary part to obtain `pv_imag`, the real part for `pv_real`).

`AbstractQAtlas` owns the generic function only; the numerical transform is the functional
sibling's job (#14 / #19), which adds a method for its own [`AbstractResponse`](@ref)
representation.  This fallback errors informatively.
"""
function principal_value_hilbert(response::AbstractResponse, ω)
    return error(
        "principal_value_hilbert not implemented for $(typeof(response)). " *
        "The principal-value Hilbert transform is a functional-sibling operation " *
        "(#14 / #19) — the numerical package must define " *
        "`principal_value_hilbert(::$(typeof(response)), ω)`.",
    )
end
export principal_value_hilbert

"""
    spectral_moment(response::AbstractResponse, n::Integer) -> Number

The `n`-th frequency moment `∫ ωⁿ f(ω) dω` of `response`: `n = 0` is the sum-rule /
normalization integral a [`SpectralSumRule`](@ref) / [`StaticFromDynamicalStructureFactor`](@ref)
check consumes (`∫A`, `∫S`), and `n = 1` the first moment a [`FSumRule`](@ref) check
consumes (`∫ω S`).

`AbstractQAtlas` owns the generic function only; the quadrature is the functional
sibling's job (#14 / #19).  This fallback errors informatively.
"""
function spectral_moment(response::AbstractResponse, n::Integer)
    return error(
        "spectral_moment not implemented for $(typeof(response)). " *
        "A spectral quadrature is a functional-sibling operation (#14 / #19) — the " *
        "numerical package must define `spectral_moment(::$(typeof(response)), n::Integer)`.",
    )
end
export spectral_moment
