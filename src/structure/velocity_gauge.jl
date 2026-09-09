# structure/velocity_gauge.jl — how a vector potential enters a lattice Hamiltonian, and what
# the current then IS. Two conventions, stated once.
#
# Neither is a fact about any model. The Peierls phase on a bond is `A · d` with `d` the hop
# distance, so it depends on the LENGTH UNIT the model measures `A` in — a choice. And the
# current is `J = -∂H/∂A` — a sign that every consumer must apply and that no static check
# can see.
#
# They live here because a calculation with a k-space arm, a real-space arm and a
# tensor-network arm needs both in all three, and three private copies is how three arms come
# to share one convention, agree with each other, and be wrong together.

"""
    PeierlsConvention(bonds_per_cell)

The length unit a vector potential is measured in, as the number of hops a unit cell spans.

`PeierlsConvention(1)` measures `A` in units of the SITE spacing, so one bond carries
`e^{iA}`. `PeierlsConvention(2)` measures it in units of a two-site CELL, so each bond spans
half a unit and carries `e^{iA/2}`.

Neither is more correct — they are different units — and the difference is invisible to most
checks: the two agree on every static quantity, on the linear response and on the
second-order response, and separate only at third order. So this is something a model
STATES; it is not something a test discovers.
"""
struct PeierlsConvention
    bonds_per_cell::Int
    function PeierlsConvention(n::Integer)
        n >= 1 || throw(ArgumentError("PeierlsConvention: need ≥1 bond per cell, got $n"))
        return new(Int(n))
    end
end
export PeierlsConvention

"""
    peierls_phase(convention, A) -> Real

The phase on ONE bond: `A / bonds_per_cell`.

A hopping becomes `exp(-i·peierls_phase(c, A))` in the forward direction and its conjugate in
the reverse. Applying the same sign to both is not a gauge transformation, and shows up as an
open chain whose energy moves with `A`.
"""
peierls_phase(c::PeierlsConvention, A::Real) = A / c.bonds_per_cell
export peierls_phase

"""
    current_from_hamiltonian_derivative(dH_dA)

`J = -∂H/∂A`, for a number, an expectation value, or an operator, given the derivative.

One negation, named because a driven current oscillates about zero and the checks that watch
it mostly take `|J|`: a flipped sign passes them and inverts every odd-order response.

Use this when the derivative is already in hand — a closed form, or a finite difference of a
measured energy. When the Hamiltonian is available as a FUNCTION of `A`, prefer
[`peierls_current`](@ref), which differentiates it rather than asking a caller to.
"""
current_from_hamiltonian_derivative(dH_dA) = -dH_dA
export current_from_hamiltonian_derivative

"""
    peierls_current(H_of_A, A) -> value

`J = -∂H/∂A` evaluated at `A` by automatic differentiation of `H_of_A`.

The same shape as [`thermal_derivative`](@ref) and the same sign convention as the rest of
the response genealogy (`M = -∂F/∂h`): a quantity IS a signed derivative of a potential, and
the seam differentiates so no consumer hand-codes one. `H_of_A` returns the energy — a number
or an expectation value — for a given vector potential.

The `A` this is differentiated at is measured in the model's own [`PeierlsConvention`](@ref),
so the length unit is the caller's to fix; this seam only supplies the derivative and the
sign.

Requires an automatic-differentiation backend to be loaded; the method is provided by the
`ForwardDiff` package extension. Without it, this throws an informative error.
"""
function peierls_current(H_of_A, A)
    return error(
        "peierls_current needs an automatic-differentiation backend — " *
        "run `using ForwardDiff` to load the AbstractQAtlas AD extension.",
    )
end
export peierls_current
