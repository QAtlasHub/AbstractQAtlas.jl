# The ForwardDiff extension: thermal_derivative evaluates the response
# genealogy by AD, cross-checked against the closed-form derivatives AND
# the supplied-derivative relations they feed.

using AbstractQAtlas
using AbstractQAtlas: thermal_derivative, check
using ForwardDiff   # activates AbstractQAtlasForwardDiffExt

@testset "single-spin free energy F(h) = −ln(2cosh βh)/β" begin
    β = 1.3
    F(h) = -log(2 * cosh(β * h)) / β
    h = 0.4
    # M = −∂F/∂h = tanh(βh)
    M = thermal_derivative(Magnetization(:z), F, h)
    @test M ≈ tanh(β * h) atol = 1e-12
    # …and it satisfies the definitional relation with the AD derivative
    @test check(
        MagnetizationResponse(); M=M, dF_dh=ForwardDiff.derivative(F, h), atol=1e-13
    )

    # χ = −∂²F/∂h² = β sech²(βh); consistent with the FDT χ = β Var(M),
    # Var(M) = 1 − tanh²(βh) = sech²(βh)
    χ = thermal_derivative(Susceptibility(:z, :z), F, h)
    @test χ ≈ β / cosh(β * h)^2 atol = 1e-10
    @test check(SusceptibilityFDT(); χ=χ, var_M=1 - tanh(β * h)^2, β=β, N=1, atol=1e-10)

    # nonlinear χ⁽²⁾ = ∂²M/∂h² = −2β² sech²(βh) tanh(βh) (the 3rd F-derivative)
    χ2 = thermal_derivative(Susceptibility(:z, :z, :z), F, h)
    @test χ2 ≈ -2β^2 * tanh(β * h) / cosh(β * h)^2 atol = 1e-9
end

@testset "tensor susceptibility: off-diagonal guard + multi-field mixed partial" begin
    # a single-field F(h) fixes only the DIAGONAL χ⁽ⁿ⁾; asking for an
    # off-diagonal component must error, not silently return the diagonal.
    @test_throws ErrorException thermal_derivative(
        Susceptibility(:x, :y, :z), h -> h^4, 0.5
    )
    # the diagonal single-field path still works: −F‴(h) = −24h
    @test thermal_derivative(Susceptibility(:x, :x, :x), h -> h^4, 0.5) ≈ -12.0 atol = 1e-10

    # multi-field: χ⁽ⁿ⁾_{α;β₁…βₙ} = −∂ⁿ⁺¹F/∂h_α∂h_{β₁}…∂h_{βₙ} over the
    # directions in indices(χ).  G = h_x h_y h_z ⇒ ∂³/∂x∂y∂z = 1 ⇒ χ = −1,
    # a value the single-field (diagonal) path can never produce.
    G(v) = v[1] * v[2] * v[3]
    @test thermal_derivative(Susceptibility(:x, :y, :z), G, [0.3, 0.7, 0.2], (:x, :y, :z)) ≈
        -1.0 atol = 1e-10
    # order-3 off-diagonal χ⁽³⁾_{w;xyz} = −∂⁴(h_w h_x h_y h_z) = −1
    G4(v) = v[1] * v[2] * v[3] * v[4]
    @test thermal_derivative(
        Susceptibility(:w, :x, :y, :z), G4, [0.1, 0.2, 0.3, 0.4], (:w, :x, :y, :z)
    ) ≈ -1.0 atol = 1e-10

    # CONSISTENCY: on the diagonal the multi-field and single-field APIs agree.
    H(v) = v[1]^4
    @test thermal_derivative(Susceptibility(:x, :x, :x), H, [0.5, 0.0, 0.0], (:x, :y, :z)) ≈
        thermal_derivative(Susceptibility(:x, :x, :x), h -> h^4, 0.5) atol = 1e-10

    # an unknown field direction errors
    @test_throws ErrorException thermal_derivative(
        Susceptibility(:x, :q), G, [0.1, 0.2, 0.3], (:x, :y, :z)
    )
end

@testset "temperature branch: S, C, U by AD" begin
    # entropy from a free-energy(T): S = −∂F/∂T
    Φ(T) = -T * log(2 * cosh(1.0 / T))
    T = 2.0
    S = thermal_derivative(ThermalEntropy(), Φ, T)
    @test check(EntropyResponse(); S=S, dF_dT=ForwardDiff.derivative(Φ, T), atol=1e-12)

    # specific heat from energy(T): C = ∂U/∂T
    U(T) = 1.5 * T                        # linear ⇒ C = 3/2
    @test thermal_derivative(SpecificHeat(), U, 2.0) ≈ 1.5 atol = 1e-12

    # Gibbs–Helmholtz: U = ∂(βF)/∂β from the βF(β) function
    βF(β) = -log(1 + exp(-β))             # two-level (ε=1): βF = −ln Z
    β0 = 0.7
    U_ad = thermal_derivative(Energy(), βF, β0)
    @test U_ad ≈ exp(-β0) / (1 + exp(-β0)) atol = 1e-12    # ⟨E⟩ of the two-level system
end

@testset "unsupported quantity falls through to the informative stub" begin
    # FreeEnergy has no AD genealogy edge ⇒ the generic stub still errors
    err = try
        thermal_derivative(FreeEnergy(), sin, 0.3)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("automatic-differentiation", err.msg)
end
