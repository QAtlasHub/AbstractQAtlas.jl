# relations/finite_size.jl: auto-discovery over the SIZES present in a bag.
#
# The size twin of `region_entropy.jl`. Keyed on a `SizeSupport`, a finite-size
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

# `Ω ~ L^-z`, so the secant of `ln Ω` against `ln L` is `-z` on a clean power law,
# for any amplitude: it is an additive constant in log space and cancels. Only
# defined for positive values, and a gap that closes exactly at one size is an
# ordinary event, so that returns no row rather than a residual of -Inf.
function _conventional_secant((L1, O1), (L2, O2))
    (O1 > 0 && O2 > 0) || return nothing
    return (log(O2) - log(O1)) / (log(L2) - log(L1))
end

# `ln(1/O) ~ L^ψ`, so it is `ln[-ln O]` that is linear in `ln L`, again for any
# amplitude. Only defined where the inner log is, which is what `0 < O < 1` says.
function _activated_secant((L1, O1), (L2, O2))
    all(O -> 0 < O < 1, (O1, O2)) || return nothing
    return (log(-log(O2)) - log(-log(O1))) / (log(L2) - log(L1))
end

"""
    finite_size_scaling_report(b::Bag; atol = 0, rtol = 0, outputlevel = 0)
        -> Vector{FiniteSizeRow}

Read every finite-size sweep in `b` against the scaling law the exponents in `b`
claim for it.

A sweep is the entries of one quantity keyed by [`at_size`](@ref). Consecutive
sizes give a secant, and that is the derivative the supplied-derivative relations
take: `d lnΩ / d lnL` against [`ConventionalFiniteSizeEnergy`](@ref) when a
[`DynamicalExponent`](@ref) is in the bag, and `d ln[-ln O] / d lnL` against
[`ActivatedFiniteSizeScaling`](@ref) when an [`ActivatedExponent`](@ref) is. Both
when both are, because that is the comparison a size sweep is usually run to
settle: the same numbers cannot obey a power of `L` and a power of `ln L`, and the
report says which they do. The secant is exact for any amplitude, which is an
additive constant in log space and cancels in a difference.

`pass` at the default `atol = 0` asks for an exact hit and measured data will not
give one; the number to read there is `residual`. `rtol` is the practical knob,
taken against the exponent being checked, so `rtol = 0.02` accepts a secant within
two percent of `ψ` or `z`.

The activated law is a statement about the TYPICAL value, so it is read only off
`Typical`-keyed data. A `DisorderAveraged` or unreduced sweep gets no activated
row: at an infinite-randomness fixed point the average is set by the rare regions
([IgloiMonthus2005](@cite), §2.2) and carries a power of `L` rather than of `ln L`,
so a secant off it would return a number that is not ψ. Its conventional row is
still reported, and so are every other sweep's rows, which is why this is a skip
and not a refusal: a bag holding typical and average data side by side is the
normal shape of an infinite-randomness study, and one sweep must not blank the
report. Pass `outputlevel = 1` to be told what was skipped and why.

```julia
b = bag(at_size(Typical(MassGap()), 16) => exp(-1.0 * sqrt(16)),
        at_size(Typical(MassGap()), 64) => exp(-1.0 * sqrt(64)),
        ActivatedExponent => 0.5)
finite_size_scaling_report(b; rtol=1e-3)
```
"""
function finite_size_scaling_report(b::Bag; atol=0, rtol=0, outputlevel=0)
    out = FiniteSizeRow[]
    sweeps = _size_sweeps(b)
    ψ = get(b, VariableKey(ActivatedExponent), nothing)
    z = get(b, VariableKey(DynamicalExponent), nothing)
    if isempty(sweeps) || (ψ === nothing && z === nothing)
        outputlevel > 0 &&
            @info "finite_size_scaling_report: nothing to read" sweeps = length(sweeps) has_ψ =
                ψ !== nothing has_z = z !== nothing
        return out
    end
    for Q in sort!(collect(keys(sweeps)); by=string)
        pts = sweeps[Q]
        if length(pts) < 2
            outputlevel > 0 &&
                @info "finite_size_scaling_report: $Q has one size, so no secant"
            continue
        end
        activated = ψ === nothing ? false : _activated_admits(Q, outputlevel)
        for (p1, p2) in zip(pts, pts[2:end])
            p1[1] == p2[1] && continue
            sizes = (p1[1], p2[1])
            if z !== nothing
                sc = _conventional_secant(p1, p2)
                if sc === nothing
                    outputlevel > 0 && @info "finite_size_scaling_report: $Q at $sizes " *
                        "is not positive, so ln Ω has nowhere to stand"
                else
                    _finite_size_scaling_row!(
                        out,
                        ConventionalFiniteSizeEnergy(),
                        Q,
                        sizes,
                        (; dlogΩ_dlogL=sc, z),
                        max(atol, rtol * abs(z)),
                    )
                end
            end
            if activated
                sa = _activated_secant(p1, p2)
                if sa === nothing
                    outputlevel > 0 && @info "finite_size_scaling_report: $Q at $sizes " *
                        "is outside 0 < O < 1, so ln[-ln O] has nowhere to stand"
                else
                    _finite_size_scaling_row!(
                        out,
                        ActivatedFiniteSizeScaling(),
                        Q,
                        sizes,
                        (; dloglogO_dlogL=sa, ψ),
                        max(atol, rtol * abs(ψ)),
                    )
                end
            end
        end
    end
    return out
end
export finite_size_scaling_report

# `ψ` in the bag makes the activated law available, and that law is about the
# typical value alone. A skip rather than an error, so one sweep cannot cost the
# rest of the bag its rows.
function _activated_admits(@nospecialize(Q::Type), outputlevel)
    Q <: Typical && return true
    if outputlevel > 0
        why = if Q <: DisorderAveraged
            "the average is set by the rare regions and carries a power of L, not of ln L"
        else
            "it says no reduction, and typical and average scale differently here"
        end
        @info "finite_size_scaling_report: no activated row for $Q, $why"
    end
    return false
end

function _finite_size_scaling_row!(out, rel, Q, sizes, vars, atol)
    r = residual(rel; vars...)
    push!(out, FiniteSizeRow(rel, Q, sizes, r, _passes(rel, r, atol)))
    return out
end
