# structure/velocity_gauge.jl — how a vector potential enters a lattice Hamiltonian, and what
# the current then IS. Definitions, stated once.
#
# The Peierls phase on the bond from site `j` to site `i` is `A · (r_i - r_j)`. Two things in
# that are conventions rather than facts, and both are invisible to every static check and to
# the linear and second-order responses:
#
#   * the LENGTH UNIT the displacement is measured in, and
#   * the sign in `J = -∂H/∂A`.
#
# They live here because a calculation with a k-space arm, a real-space arm and a
# tensor-network arm needs both in all three, and three private copies is how three arms come
# to share one convention, agree with each other, and be wrong together.
#
# `A` is a VECTOR. The dimension is carried in the type rather than left to a bare `Real`, so
# a seam that only handles one dimension has to say so instead of accepting a number and
# meaning something narrower than it looks.
#
# WHERE THIS SITS IN THE EXISTING LAYERS. Three of them already carry a piece of this and
# none of them carried the phase:
#
#   * `SpatialDirection` (core/indices.jl) names the index `μ` that currents and transport
#     tensors carry — `σ_μν`, `DynamicalConductivity(:x, :y)`. `VectorPotential{N}` is a
#     VALUE in that same space, so `N` counts `SpatialDirection` slots. The existing rule
#     there is that "its range is set by the model, so the interface names the space without
#     enumerating it"; `N` is that range, made concrete once a caller supplies a field.
#   * `Region` (core/region.jl) is the SUPPORT a quantity is evaluated on, and is
#     deliberately dimension-agnostic: a site is any hashable label. That is the right shape
#     for WHICH sites, and it is not the shape for the displacement BETWEEN two of them —
#     `A · (r_i - r_j)` needs positions in a vector space, which a label set does not have.
#     The two are complementary: `Region` says where, `VectorPotential` says along what.
#   * `DynamicalConductivity` and `CurrentCorrelation` (core/quantities.jl) already index by
#     `SpatialDirection`. `J = -∂H/∂A` is the operator whose correlations those quantities
#     ARE, so the sign fixed here is the sign their Kubo expressions inherit.
#
# The layer added here is the one between them: the field, its contraction with a
# displacement, and the sign of the current that contraction defines.

"""
    VectorPotential(components::Real...)
    VectorPotential(a::Real)                # one dimension

The optical vector potential `A`, with its dimension in the type.

`N` IS THE DIMENSION OF THE SPACE DISPLACEMENTS LIVE IN — the ambient space the sites are
embedded in — because [`peierls_phase`](@ref) contracts `A` with `r_i - r_j`. It is not a
property of the site SET, and in particular it is not a Hausdorff dimension: a Sierpiński
gasket drawn in the plane has sites whose set is `log3/log2`-dimensional and displacements
that are ordinary 2-vectors, so its `A` is `VectorPotential{2}`.

For a cut-and-project quasicrystal it is the PARALLEL dimension (`D_par` in
QuasiCrystal.jl's `cut_and_project_dimensions`), not the hyperspace `D_hyper`: a Fibonacci
chain has `D_par = 1`, `D_hyper = 2`, and its optical `A` is `VectorPotential{1}` because the
hops it multiplies are displacements in the physical line. A drive along the PERPENDICULAR
directions is a phason, not a vector potential, and does not belong in this type — conflating
the two is the mistake `N` exists to make impossible to write.

A bare number could say none of this.
"""
struct VectorPotential{N,T<:Real}
    components::NTuple{N,T}
end
VectorPotential(components::Real...) = VectorPotential(promote(components...))
VectorPotential(a::Real) = VectorPotential((a,))
export VectorPotential

"""
    dimension(A::VectorPotential) -> Int

How many spatial components `A` has.
"""
dimension(::VectorPotential{N}) where {N} = N
export dimension

Base.getindex(A::VectorPotential, i::Int) = A.components[i]
Base.length(A::VectorPotential) = length(A.components)
Base.iterate(A::VectorPotential, s...) = iterate(A.components, s...)

"""
    peierls_phase(A::VectorPotential{N}, displacement::NTuple{N,<:Real}) -> Real

`A · d`, the phase a hopping picks up across a bond whose endpoints differ by `displacement`.

The hopping becomes `exp(-i·peierls_phase(A, d))` in the forward direction and its conjugate
in the reverse. Applying the same sign to both is not a gauge transformation, and shows up as
an open chain whose energy moves with `A`.

`displacement` is measured in whatever length unit the model uses — see
[`PeierlsConvention`](@ref) for the one-dimensional case, where every bond spans the same
distance and the unit is the whole of the choice.
"""
function peierls_phase(A::VectorPotential{N}, displacement::NTuple{N,<:Real}) where {N}
    return sum(A.components .* displacement)
end
export peierls_phase

"""
    PeierlsConvention(bonds_per_cell)

The one-dimensional shortcut: a uniform chain whose every bond spans `1 / bonds_per_cell`, so
naming that integer names the length unit.

`PeierlsConvention(1)` measures `A` in units of the SITE spacing, so one bond carries
`e^{iA}`. `PeierlsConvention(2)` measures it in units of a two-site CELL, so each bond spans
half a unit and carries `e^{iA/2}`.

Neither is more correct — they are different units — and the difference separates only at
third order in the drive. So this is something a model STATES; it is not something a test
discovers.

Beyond one dimension the displacements are not all equal and there is no single integer to
name; pass them to [`peierls_phase`](@ref) directly.
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
    bond_displacement(c::PeierlsConvention) -> Real

`1 / bonds_per_cell`, the displacement every bond of the uniform chain spans.

Separated from [`peierls_phase`](@ref) because it is the model-dependent half: the phase is
always `A · d`, and this is what `d` happens to be in one dimension.
"""
bond_displacement(c::PeierlsConvention) = 1 / c.bonds_per_cell
export bond_displacement

"""
    peierls_phase(c::PeierlsConvention, A::VectorPotential{1}) -> Real

The one-dimensional phase, `A / bonds_per_cell`. The general `A · d` with the displacement
[`bond_displacement`](@ref) supplies.
"""
function peierls_phase(c::PeierlsConvention, A::VectorPotential{1})
    return peierls_phase(A, (bond_displacement(c),))
end

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
    peierls_current(H_of_A, A::VectorPotential) -> value

`J = -∂H/∂A` at `A`, by automatic differentiation of `H_of_A`.

The same shape as [`thermal_derivative`](@ref) and the same sign convention as the rest of
the response genealogy (`M = -∂F/∂h`): a quantity IS a signed derivative of a potential, and
the seam differentiates so no consumer hand-codes one. `H_of_A` takes a
[`VectorPotential`](@ref) and returns the energy.

Currently implemented for `VectorPotential{1}` only, by forward mode; a higher dimension
needs a gradient and is refused rather than silently reduced to one component. The length
unit `A` is measured in stays the model's business — this seam supplies the derivative and
the sign, nothing else.

Requires an automatic-differentiation backend; the method is provided by the `ForwardDiff`
package extension. Without it, this throws an informative error.
"""
function peierls_current(H_of_A, A::VectorPotential{N}) where {N}
    N == 1 && return error(
        "peierls_current needs an automatic-differentiation backend — " *
        "run `using ForwardDiff` to load the AbstractQAtlas AD extension.",
    )
    return error(
        "peierls_current: no method for a $(N)-dimensional vector potential. The " *
        "one-dimensional case differentiates in forward mode; a gradient needs a " *
        "reverse-mode route, which is not written yet. Reducing to one component here " *
        "would silently answer a different question.",
    )
end
export peierls_current
