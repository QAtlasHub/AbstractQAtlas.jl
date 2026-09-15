# src/derivative_routes.jl: HOW a genealogy edge is evaluated.
#
# structure/response.jl says WHICH derivative a response function is
# (`derivative_edge`, `derivative_order`). Until now there was exactly one way
# to evaluate it, nested forward-mode AD in a package extension, so a project
# holding a measured potential rather than a differentiable one could not
# traverse the genealogy at all.
#
# The genealogy work is the same whichever way the derivative is taken, so it
# lives here once, in `_genealogy_derivative`, parameterised by an `nth`
# operator. The ForwardDiff extension supplies the AD `nth` and keeps its own
# entry points; the finite-difference routes below need no backend.
#
# Finite differences were already in this package three times, written inline
# and unnamed: relations/region_entropy.jl's entropy secant and the two
# log-log secants in relations/finite_size.jl. Those stay where they are (they
# difference two MEASURED points, not a function), but they are why a route is
# a declared thing here rather than an argument spelled out at each call.

"""
    DerivativeRoute

How a derivative is evaluated: [`AutoDiff`](@ref), [`CentralDifference`](@ref),
[`Richardson`](@ref).

Open in the same way [`Convention`](@ref) is: a new route is a subtype plus an
[`nth_derivative`](@ref) method, and the genealogy layer above it does not
change. Which route a project uses is a property of what it measured, not of
the physics, which is why it is named and passed rather than fixed here.
"""
abstract type DerivativeRoute end
export DerivativeRoute

"""
    AutoDiff() <: DerivativeRoute

Nested forward-mode automatic differentiation: exact to machine precision, and
the only route that needs the potential to be a differentiable Julia function.
Requires the ForwardDiff extension to be loaded.
"""
struct AutoDiff <: DerivativeRoute end
export AutoDiff

"""
    CentralDifference(h::Real) <: DerivativeRoute

The symmetric difference quotient at step `h`, applied `n` times for an `n`-th
derivative. Second-order accurate, and needs only that the potential can be
called at points near `x`.

`h` has no default. The error is `O(h²) + O(ε/hⁿ)`, so the best step depends on
the order taken and on how noisy the potential is, and a default here would be
a number chosen without either. [`observed_order`](@ref) is how to tell whether
the `h` in hand is on the truncation side or the roundoff side.
"""
struct CentralDifference <: DerivativeRoute
    h::Float64
    function CentralDifference(h::Real)
        h > 0 || throw(ArgumentError("CentralDifference: h must be positive; got $h"))
        return new(Float64(h))
    end
end
export CentralDifference

"""
    Richardson(h::Real; levels::Int = 3) <: DerivativeRoute

[`CentralDifference`](@ref) evaluated at `h, h/2, …` and extrapolated, removing
`levels - 1` orders of the truncation error.

The accurate route when there is no AD backend: on a smooth potential it
reaches AD to several digits where a single central difference at the same `h`
reaches two or three.
"""
struct Richardson <: DerivativeRoute
    h::Float64
    levels::Int
    function Richardson(h::Real; levels::Int=3)
        h > 0 || throw(ArgumentError("Richardson: h must be positive; got $h"))
        levels >= 2 ||
            throw(ArgumentError("Richardson: levels must be at least 2; got $levels"))
        return new(Float64(h), levels)
    end
end
export Richardson

"""
    nth_derivative(route::DerivativeRoute, f, x, n::Integer) -> value

The `n`-th derivative of the scalar function `f` at `x`, taken along `route`.
`n = 0` is `f(x)` on every route.

This is the one method a new route has to define.
"""
function nth_derivative(route::DerivativeRoute, f, x, n::Integer)
    return error("nth_derivative: no method for $(typeof(route)).")
end
export nth_derivative

function nth_derivative(::AutoDiff, f, x, n::Integer)
    return error(
        "nth_derivative(AutoDiff(), ...) needs an automatic-differentiation " *
        "backend: run `using ForwardDiff` to load the AbstractQAtlas AD extension, " *
        "or take a finite-difference route (CentralDifference / Richardson), which " *
        "needs none.",
    )
end

_central(f, x, h) = (f(x + h) - f(x - h)) / (2h)

function nth_derivative(route::CentralDifference, f, x, n::Integer)
    n >= 0 || throw(ArgumentError("nth_derivative: n must be non-negative; got $n"))
    n == 0 && return f(x)
    return _central(y -> nth_derivative(route, f, y, n - 1), x, route.h)
end

function nth_derivative(route::Richardson, f, x, n::Integer)
    n >= 0 || throw(ArgumentError("nth_derivative: n must be non-negative; got $n"))
    n == 0 && return f(x)
    # Neville on the halved steps: each column removes one more order of h².
    t = [
        nth_derivative(CentralDifference(route.h / 2^k), f, x, n) for
        k in 0:(route.levels - 1)
    ]
    for j in 1:(route.levels - 1)
        w = 4.0^j
        t = [(w * t[i + 1] - t[i]) / (w - 1) for i in 1:(length(t) - 1)]
    end
    return only(t)
end

"""
    observed_order(route::DerivativeRoute, f, x, n::Integer) -> Float64

The convergence order the route actually shows on `f` at `x`, from the values at
`h`, `h/2` and `h/4`: `log2(|D(h) - D(h/2)| / |D(h/2) - D(h/4)|)`.

The number to look at before trusting a step, rather than a tolerance guessed in
advance. A central difference on a smooth potential returns close to 2; a value
well below that means `h` has reached the roundoff side, and a value near 0 means
`f` is not smooth at `x`. Returns `NaN` when the two differences are both zero,
which is the step being so small that the quotient stopped moving.

Meaningful only while the successive differences are above roundoff. A route that
has already reached machine precision, which [`Richardson`](@ref) does on a smooth
potential, is differencing noise and reports a number with no order in it.

Defined for the step-carrying routes; [`AutoDiff`](@ref) has no step to halve.
"""
function observed_order(route::DerivativeRoute, f, x, n::Integer)
    return error("observed_order: $(typeof(route)) carries no step to halve.")
end
export observed_order

_with_step(r::CentralDifference, h) = CentralDifference(h)
_with_step(r::Richardson, h) = Richardson(h; levels=r.levels)
_step(r::CentralDifference) = r.h
_step(r::Richardson) = r.h

function observed_order(route::Union{CentralDifference,Richardson}, f, x, n::Integer)
    h = _step(route)
    d = [nth_derivative(_with_step(route, h / 2^k), f, x, n) for k in 0:2]
    a, b = abs(d[2] - d[1]), abs(d[3] - d[2])
    (a == 0 && b == 0) && return NaN
    b == 0 && return Inf
    return log2(a / b)
end

# ── the genealogy, written once ──────────────────────────────────────────
#
# Shared with the ForwardDiff extension, which passes its own `nth`. The net
# sign is -1 across this tree: the `F -> response` edge is `-∂/∂field` and the
# intra-response edges are `+∂/∂field`.
function _genealogy_derivative(q::AbstractQuantity, F, x, nth)
    e = derivative_edge(q)
    e === nothing && error(
        "thermal_derivative: $(typeof(q)) is not a response function " *
        "(no derivative_edge), so there is nothing to differentiate.",
    )
    root = potential_root(q)
    (root === FreeEnergy || root === GrandPotential) || error(
        "thermal_derivative: the generic genealogy path handles the standard response " *
        "potentials (FreeEnergy `M=-∂F/∂h`, GrandPotential `N=-∂Ω/∂μ`); $(typeof(q)) " *
        "roots at $root, so use its explicit method (e.g. Energy/SpecificHeat).",
    )
    return -nth(F, x, derivative_order(q, e.field()))
end

"""
    thermal_derivative(quantity, potential, x, route::DerivativeRoute) -> value

[`thermal_derivative`](@ref) along an explicit [`DerivativeRoute`](@ref), so a
project that measured its potential on a grid can traverse the response
genealogy without an AD backend.

`route = AutoDiff()` is the three-argument method, and needs the extension. The
finite-difference routes need nothing.

```julia
F(h) = -log(2cosh(h))
thermal_derivative(Magnetization(:z), F, 0.3, Richardson(1e-2))   # tanh(0.3)
```
"""
function thermal_derivative(q::AbstractQuantity, F, x::Number, route::DerivativeRoute)
    route isa AutoDiff && return thermal_derivative(q, F, x)
    return _genealogy_derivative(q, F, x, (g, y, n) -> nth_derivative(route, g, y, n))
end

# The two members whose potential is not the root: `C = ∂U/∂T` is reached through
# `U`, and `U = ∂(βF)/∂β` takes `βF`. Both are a plain first derivative of the
# function passed, with no sign flip, which is why they cannot go through the
# generic path above.
function thermal_derivative(::SpecificHeat, U, T::Number, route::DerivativeRoute)
    route isa AutoDiff && return thermal_derivative(SpecificHeat(), U, T)
    return nth_derivative(route, U, T, 1)
end
function thermal_derivative(::Energy, βF, β::Number, route::DerivativeRoute)
    route isa AutoDiff && return thermal_derivative(Energy(), βF, β)
    return nth_derivative(route, βF, β, 1)
end

# A single-field potential fixes only the DIAGONAL susceptibility; the same guard
# the AD path carries, so a route change cannot turn a refusal into a wrong number.
function thermal_derivative(χ::Susceptibility, F, h::Number, route::DerivativeRoute)
    route isa AutoDiff && return thermal_derivative(χ, F, h)
    idx = indices(χ)
    all(==(idx[1]), idx) || error(
        "thermal_derivative: with a single-field function only the DIAGONAL χ⁽ⁿ⁾ " *
        "(all indices equal) is defined; got off-diagonal $(idx). Pass a multi-field " *
        "potential F(h⃗) and the field-component ordering.",
    )
    return -nth_derivative(route, F, h, response_order(χ) + 1)
end

# The `n` the route is asked for, so a report on a third-order response halves its
# step against the third derivative and not the first.
function _route_order(q::AbstractQuantity)
    e = derivative_edge(q)
    e === nothing && return 1
    return derivative_order(q, e.field())
end

"""
    DerivativeRouteRow

One row of a [`derivative_report`](@ref): the `route`, the `value` it returned,
and the `order` it showed ([`observed_order`](@ref); `NaN` for a route with no
step to halve).
"""
struct DerivativeRouteRow
    route::DerivativeRoute
    value::Float64
    order::Float64
end
export DerivativeRouteRow

"""
    derivative_report(quantity, potential, x, routes) -> Vector{DerivativeRouteRow}

`quantity` evaluated along each of `routes`, so the routes can be compared
rather than trusted one at a time.

Two routes agreeing is evidence a single route cannot give: a step on the
roundoff side and a step on the truncation side both return a number, and only
the spread between them says which. A route that throws is reported as `NaN`
rather than aborting the sweep, since the usual reason is a missing backend and
the other rows are still the answer.
"""
function derivative_report(q::AbstractQuantity, F, x::Number, routes)
    out = DerivativeRouteRow[]
    for r in routes
        v = try
            Float64(thermal_derivative(q, F, x, r))
        catch
            NaN
        end
        o = try
            Float64(observed_order(r, F, x, _route_order(q)))
        catch
            NaN
        end
        push!(out, DerivativeRouteRow(r, v, o))
    end
    return out
end
export derivative_report
