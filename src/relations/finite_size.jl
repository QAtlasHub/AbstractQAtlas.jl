# relations/finite_size.jl — auto-discovery over the SIZES present in a bag.
#
# The size twin of `region_entropy.jl`: keyed on a `SizeSupport`, a finite-size
# sweep becomes something a bag can hold and this layer can read, so the
# supplied-derivative relations get their derivative from the data rather than
# from the caller. Two sizes are enough, and on a clean power law the secant is
# the derivative exactly.

"""
    FiniteSizeRow

One row of a [`finite_size_scaling_report`](@ref): the `relation` matched, the
`quantity` whose sweep it was read from, the `sizes` the secant was taken
between, its [`residual`](@ref), and `pass`.
"""
struct FiniteSizeRow
    relation::AbstractRelation
    quantity::Type
    sizes::Tuple{Any,Any}
    residual::Number
    pass::Bool
end
export FiniteSizeRow

# The sweeps in a bag, grouped by quantity type and sorted by size.
function _size_sweeps(b::Bag)
    out = Dict{Type,Vector{Tuple{Any,Any}}}()
    for (k, v) in b
        k.support isa SizeSupport || continue
        push!(get!(out, k.type, Tuple{Any,Any}[]), (k.support.size, v))
    end
    for v in values(out)
        sort!(v; by=first)
    end
    return out
end

# `Ω ~ L^-z`, so the secant of `ln Ω` against `ln L` is `-z` on a clean power law.
_conventional_secant((L1, O1), (L2, O2)) = (log(O2) - log(O1)) / (log(L2) - log(L1))

# `ln(1/O) ~ L^ψ`, so it is `ln[-ln O]` that is linear in `ln L`. Only defined
# where the inner log is, which is what `0 < O < 1` says.
function _activated_secant((L1, O1), (L2, O2))
    all(O -> 0 < O < 1, (O1, O2)) || return nothing
    return (log(-log(O2)) - log(-log(O1))) / (log(L2) - log(L1))
end

"""
    finite_size_scaling_report(b::Bag; atol = 1e-8) -> Vector{FiniteSizeRow}

Read every finite-size sweep in `b` against the scaling law the exponents in `b`
claim for it.

A sweep is the entries of one quantity keyed by [`at_size`](@ref). Consecutive
sizes give a secant, and that is the derivative the supplied-derivative
relations take: `d lnΩ / d lnL` against [`ConventionalFiniteSizeEnergy`](@ref)
when a [`DynamicalExponent`](@ref) is in the bag, and `d ln[-ln O] / d lnL`
against [`ActivatedFiniteSizeScaling`](@ref) when an [`ActivatedExponent`](@ref)
is. Both, when both are, because that is the comparison a finite-size sweep is
usually run to settle: the same numbers cannot obey a power of `L` and a power
of `ln L` at once, and the report says which one they do.

The activated law is a statement about the TYPICAL value, so it is read only off
`Typical`-keyed data. A sweep of `DisorderAveraged` or of a bare quantity is
refused rather than skipped while `ψ` is present: the average of a gap at an
infinite-randomness fixed point is set by the rare regions and carries a
different exponent, so fitting it to `L^ψ` returns a number that looks like an
answer. Refusing is the same call [`singular_form`](@ref) makes.

```julia
b = bag(at_size(Typical(MassGap()), 16) => exp(-1.0 * sqrt(16)),
        at_size(Typical(MassGap()), 64) => exp(-1.0 * sqrt(64)),
        ActivatedExponent => 0.5)
finite_size_scaling_report(b)
```
"""
function finite_size_scaling_report(b::Bag; atol=1e-8)
    out = FiniteSizeRow[]
    sweeps = _size_sweeps(b)
    isempty(sweeps) && return out
    ψ = get(b, VariableKey(ActivatedExponent), nothing)
    z = get(b, VariableKey(DynamicalExponent), nothing)
    (ψ === nothing && z === nothing) && return out
    for Q in sort!(collect(keys(sweeps)); by=string)
        pts = sweeps[Q]
        length(pts) >= 2 || continue
        ψ === nothing || _refuse_untyped_sweep(Q)
        for (p1, p2) in zip(pts, pts[2:end])
            p1[1] == p2[1] && continue
            sizes = (p1[1], p2[1])
            if z !== nothing
                _finite_size_scaling_row!(
                    out,
                    ConventionalFiniteSizeEnergy(),
                    Q,
                    sizes,
                    (; dlogΩ_dlogL=_conventional_secant(p1, p2), z),
                    atol,
                )
            end
            if ψ !== nothing
                s = _activated_secant(p1, p2)
                s === nothing || _finite_size_scaling_row!(
                    out,
                    ActivatedFiniteSizeScaling(),
                    Q,
                    sizes,
                    (; dloglogO_dlogL=s, ψ),
                    atol,
                )
            end
        end
    end
    return out
end
export finite_size_scaling_report

# `ψ` in the bag makes the activated law available, and the activated law is
# about the typical value alone.
function _refuse_untyped_sweep(@nospecialize(Q::Type))
    Q <: Typical && return nothing
    Q <: DisorderAveraged && error(
        "finite_size_scaling_report: the sweep is DisorderAveraged and the bag " *
        "carries an ActivatedExponent. At an infinite-randomness fixed point the " *
        "average is set by the rare regions and does not follow L^ψ, so a secant " *
        "off it would return a number that is not ψ. Key the sweep as Typical(...), " *
        "or drop ActivatedExponent to read it conventionally.",
    )
    return error(
        "finite_size_scaling_report: the sweep is $Q, with no reduction, and the " *
        "bag carries an ActivatedExponent. Typical and average scale differently " *
        "there, so say which this is: Typical(...) or DisorderAveraged(...).",
    )
end

function _finite_size_scaling_row!(out, rel, Q, sizes, vars, atol)
    r = residual(rel; vars...)
    push!(out, FiniteSizeRow(rel, Q, sizes, r, isapprox(r, zero(r); atol=atol)))
    return out
end
