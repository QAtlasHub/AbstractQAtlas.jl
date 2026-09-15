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
    backend_package(route) === nothing || throw(MissingRouteBackend(route))
    return error("nth_derivative: no method for $(typeof(route)).")
end
export nth_derivative

"""
    MissingRouteBackend(route) <: Exception

Thrown when a [`DerivativeRoute`](@ref) needs a package extension that is not
loaded.

Its own type, because [`derivative_report`](@ref) has to tell it from every other
way a route can fail. Catching `Exception` there would turn a diagnosed refusal,
an off-diagonal susceptibility or a potential evaluated outside its domain, into
the same `NaN` row as an unloaded backend.
"""
struct MissingRouteBackend <: Exception
    route::DerivativeRoute
end
export MissingRouteBackend

function Base.showerror(io::IO, e::MissingRouteBackend)
    return print(
        io,
        "MissingRouteBackend: $(typeof(e.route)) needs the $(backend_package(e.route)) ",
        "extension, which is not loaded. Run `using $(backend_package(e.route))`, or ",
        "take a finite-difference route (CentralDifference / Richardson), which needs ",
        "no backend.",
    )
end

"""
    backend_package(route::DerivativeRoute) -> Union{Symbol,Nothing}

The package whose extension supplies `route`'s [`nth_derivative`](@ref), or
`nothing` for a route that needs none.

A route declaring one and finding no method gets
[`MissingRouteBackend`](@ref) rather than a bare "no method", which is the
difference between "install this" and "this route does not exist". Declared here
and not in the extension: the point is to answer when the extension is ABSENT.
"""
backend_package(::DerivativeRoute) = nothing
backend_package(::AutoDiff) = :ForwardDiff
export backend_package

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
    step_size(route::DerivativeRoute) -> Union{Real,Nothing}
    with_step_size(route::DerivativeRoute, h::Real) -> DerivativeRoute

The step `route` takes, and the same route at a different step. `nothing` means
the route has no step, which is what [`AutoDiff`](@ref) reports.

Part of the route contract alongside [`nth_derivative`](@ref), and the pair
[`observed_order`](@ref) needs. A route that carries a step and defines neither
gets the honest refusal rather than the false claim that it has no step, which is
what a closed `Union` over the routes that happened to exist would have told it.
"""
step_size(::DerivativeRoute) = nothing
export step_size

function with_step_size(route::DerivativeRoute, h::Real)
    return error(
        "with_step_size: $(typeof(route)) defines no `with_step_size`. A route that " *
        "reports a `step_size` needs one, so `observed_order` can halve it.",
    )
end
export with_step_size

step_size(r::CentralDifference) = r.h
step_size(r::Richardson) = r.h
with_step_size(::CentralDifference, h::Real) = CentralDifference(h)
with_step_size(r::Richardson, h::Real) = Richardson(h; levels=r.levels)

"""
    observed_order(route::DerivativeRoute, f, x, n::Integer) -> Float64

The convergence order the route actually shows on `f` at `x`, from the values at
`h`, `h/2` and `h/4`: `log2(|D(h) - D(h/2)| / |D(h/2) - D(h/4)|)`.

The number to look at before trusting a step, rather than a tolerance guessed in
advance. A central difference on a smooth potential returns close to 2. Anything
else says the step or the potential is not what the route assumed, and the value
does not identify which: `h` on the roundoff side and a non-smooth `f` both land
off 2, and a kink gives exactly 1 rather than anything near 0.

`Inf` when only the second difference vanishes and `NaN` when both do, which is
the quotient having stopped moving between halvings.

Meaningful only while the successive differences are above roundoff. A route that
has already reached machine precision, which [`Richardson`](@ref) does on a smooth
potential, is differencing noise and reports a number with no order in it.

Defined for any route reporting a [`step_size`](@ref); [`AutoDiff`](@ref) reports
`nothing` and is refused.
"""
function observed_order(route::DerivativeRoute, f, x, n::Integer)
    h = step_size(route)
    h === nothing && error(
        "observed_order: $(typeof(route)) reports no `step_size`, so there is no step " *
        "to halve.",
    )
    d = [nth_derivative(with_step_size(route, h / 2^k), f, x, n) for k in 0:2]
    a, b = abs(d[2] - d[1]), abs(d[3] - d[2])
    (a == 0 && b == 0) && return NaN
    b == 0 && return Inf
    return log2(a / b)
end
export observed_order

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
    return _genealogy_derivative(q, F, x, (g, y, n) -> nth_derivative(route, g, y, n))
end

# The two members whose potential is not the root: `C = ∂U/∂T` is reached through
# `U`, and `U = ∂(βF)/∂β` takes `βF`. Both are a plain first derivative of the
# function passed, with no sign flip, which is why they cannot go through the
# generic path above.
function thermal_derivative(::SpecificHeat, U, T::Number, route::DerivativeRoute)
    return nth_derivative(route, U, T, 1)
end
function thermal_derivative(::Energy, βF, β::Number, route::DerivativeRoute)
    return nth_derivative(route, βF, β, 1)
end

# A single-field potential fixes only the DIAGONAL susceptibility: an off-diagonal
# component is a mixed partial in distinct field directions. Shared with the
# extension for the same reason `_genealogy_derivative` is, so a route change
# cannot turn a refusal into a wrong number.
function _susceptibility_derivative(χ::Susceptibility, F, h, nth)
    idx = indices(χ)
    all(==(idx[1]), idx) || error(
        "thermal_derivative: with a single-field function only the DIAGONAL χ⁽ⁿ⁾ " *
        "(all indices equal) is defined; got off-diagonal $(idx). Pass a multi-field " *
        "potential F(h⃗) and the field-component ordering.",
    )
    return -nth(F, h, response_order(χ) + 1)
end

function thermal_derivative(χ::Susceptibility, F, h::Number, route::DerivativeRoute)
    return _susceptibility_derivative(χ, F, h, (g, y, n) -> nth_derivative(route, g, y, n))
end

# The `n` the route is asked for, so a report on a third-order response halves its
# step against the third derivative and not the first. Callers check the edge first:
# there is no order to report for a quantity that is not a derivative of anything.
function _route_order(q::AbstractQuantity)
    e = derivative_edge(q)
    e === nothing &&
        error("_route_order: $(typeof(q)) has no derivative_edge, so it has no order.")
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
    # Refused here rather than per row: `_route_order` would fall back to 1 and the
    # order column would report a textbook 2.0 beside a value the same call refused,
    # which reads as "the numerics are fine, only the value failed".
    derivative_edge(q) === nothing && error(
        "derivative_report: $(typeof(q)) is not a response function (no " *
        "derivative_edge), so there is no derivative for a route to take.",
    )
    n = _route_order(q)
    out = DerivativeRouteRow[]
    for r in routes
        v = try
            Float64(thermal_derivative(q, F, x, r))
        catch e
            e isa MissingRouteBackend || rethrow()
            NaN
        end
        o = try
            Float64(observed_order(r, F, x, n))
        catch e
            (e isa MissingRouteBackend || e isa ErrorException) || rethrow()
            NaN
        end
        push!(out, DerivativeRouteRow(r, v, o))
    end
    return out
end
export derivative_report
