# src/autodiff.jl — automatic-differentiation entry point (stub).
#
# The supplied-derivative relations (S = −∂F/∂T, χ = ∂M/∂h, U = ∂(βF)/∂β,
# the Maxwell relations, …) take a derivative *value*.  `thermal_derivative`
# evaluates that value from the underlying potential FUNCTION by
# automatic differentiation, structured by the response genealogy
# (`derivative_edge`): the genealogy declares that a quantity is a
# derivative of a potential; this computes it.
#
# The implementation lives in a package extension so AD is an optional
# dependency — `residual`/`check`/`solve` and the forms are pure
# arithmetic and already differentiate through any AD backend without it.

"""
    thermal_derivative(quantity, potential, x) -> value
    thermal_derivative(χ::Susceptibility, F, h⃗::AbstractVector, components) -> value

The value of `quantity` as the appropriate derivative of the `potential`
function evaluated at the point `x`, via automatic differentiation — the
AD realization of the response genealogy (`derivative_edge`):

| `quantity` | `potential` | result |
|---|---|---|
| `Magnetization(α)` | `F(h)` | `M_α = −∂F/∂h` |
| `Susceptibility(α, β₁…βₙ)` | `F(h)` | `χ⁽ⁿ⁾ = −∂ⁿ⁺¹F/∂hⁿ⁺¹` (**diagonal only** — all indices equal) |
| `ThermalEntropy()` | `F(T)` | `S = −∂F/∂T` |
| `SpecificHeat()` | `U(T)` | `C = ∂U/∂T` |
| `Energy()` | `βF(β)` | `U = ∂(βF)/∂β` (Gibbs–Helmholtz) |

A single-field `F(h)` fixes only the **diagonal** susceptibility (every
index equal); an off-diagonal component is a mixed partial in distinct
field directions and errors (rather than silently returning the
diagonal).  For the full tensor component pass a **multi-field**
potential `F(h⃗)` and the field-direction ordering `components`:

`χ⁽ⁿ⁾_{α;β₁…βₙ} = −∂ⁿ⁺¹F / ∂h_α ∂h_{β₁} … ∂h_{βₙ}`

(the response index `α` is included; the diagonal reproduces the
single-field result).

Requires an automatic-differentiation backend to be loaded; the methods
are provided by the `ForwardDiff` package extension.  Without it, this
throws an informative error.

```julia
using ForwardDiff
F(h) = -log(2cosh(h)) / β                       # single-spin free energy
thermal_derivative(Magnetization(:z), F, 0.3)   # M = tanh(0.3·β)·…  (= −F'(0.3))

G(h⃗) = h⃗[1] * h⃗[2] * h⃗[3]                       # a cross-field free energy
thermal_derivative(Susceptibility(:x, :y, :z), G, [0.0, 0.0, 0.0], (:x, :y, :z))  # −1
```
"""
function thermal_derivative(quantity::AbstractQuantity, potential, x)
    return error(
        "thermal_derivative needs an automatic-differentiation backend — " *
        "run `using ForwardDiff` to load the AbstractQAtlas AD extension.",
    )
end
export thermal_derivative
