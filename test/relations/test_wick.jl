# Wick's theorem vs INDEPENDENT many-body Fock-space ED.
#
# A 3-site free-fermion chain is solved exactly in the full 8-dimensional
# Fock space (with explicit Jordan–Wigner sign bookkeeping).  Both the
# 2-point matrix G AND the 4-point functions are computed from the
# many-body thermal state directly; Wick's determinant factorization of
# the 4-points from G is then a genuinely independent identity check
# (ED trace vs determinant formula).

using AbstractQAtlas
using LinearAlgebra
using Test

# ── Fock-space machinery (3 sites, 8 states, JW convention) ─────────────
# |n⟩ = (c†₁)^{n₁}(c†₂)^{n₂}(c†₃)^{n₃}|0⟩, bit k of the integer label is
# the occupation of site k.  c_i |n⟩ carries sign (−1)^{Σ_{k<i} n_k}.
const NS = 3
const DIM = 2^NS

function annihilator(i::Int)
    c = zeros(Float64, DIM, DIM)
    for n in 0:(DIM - 1)
        if (n >> (i - 1)) & 1 == 1
            sign = iseven(count_ones(n & ((1 << (i - 1)) - 1))) ? 1.0 : -1.0
            c[n - (1 << (i - 1)) + 1, n + 1] = sign
        end
    end
    return c
end
const CS = [annihilator(i) for i in 1:NS]
const CDS = [Matrix(c') for c in CS]

@testset "Wick vs Fock-space ED (3-site tight binding, thermal)" begin
    t = 1.0
    H = zeros(DIM, DIM)
    for i in 1:(NS - 1)
        H .+= -t .* (CDS[i] * CS[i + 1] .+ CDS[i + 1] * CS[i])
    end

    for β in (0.7, 2.5)
        ρ = exp(-β .* H)
        ρ ./= tr(ρ)

        # 2-point matrix from the many-body state
        G = [tr(ρ * CDS[i] * CS[j]) for i in 1:NS, j in 1:NS]

        @testset "β = $β" begin
            # 2-point round trip through wick_contraction (n = 1 determinant)
            for i in 1:NS, j in 1:NS
                @test wick_contraction(G, [i], [j]) ≈ G[i, j] atol = 1e-12
            end

            # densities and density-density correlations: ED vs Wick
            for i in 1:NS, j in 1:NS
                nn_ed = tr(ρ * (CDS[i] * CS[i]) * (CDS[j] * CS[j]))
                @test wick_density_correlation(G, i, j) ≈ nn_ed atol = 1e-12
            end

            # a genuinely off-diagonal 4-point: ⟨c†₁ c†₂ c₂ c₃⟩
            #   operator string  c†_{cr₁} c†_{cr₂} c_{an₂} c_{an₁}
            #   with cr = [1,2], an = [3,2]  (an₂ = 2, an₁ = 3)
            val_ed = tr(ρ * CDS[1] * CDS[2] * CS[2] * CS[3])
            @test wick_contraction(G, [1, 2], [3, 2]) ≈ val_ed atol = 1e-12

            # and a 6-point: ⟨c†₁c†₂c†₃ c₃c₂c₁⟩ = ⟨n₁n₂n₃⟩
            val6_ed = tr(ρ * (CDS[1] * CS[1]) * (CDS[2] * CS[2]) * (CDS[3] * CS[3]))
            @test wick_contraction(G, [1, 2, 3], [1, 2, 3]) ≈ val6_ed atol = 1e-12
        end
    end
end

@testset "input validation" begin
    G = [0.5 0.0; 0.0 0.5]
    @test_throws ErrorException wick_contraction(G, [1, 2], [1])
end

@testset "anomalous / BdG Wick: the Pfaffian" begin
    using AbstractQAtlas: wick_pfaffian
    using LinearAlgebra: det
    # a 4×4 antisymmetric matrix: Pf = A12 A34 − A13 A24 + A14 A23 (definition)
    a12, a13, a14, a23, a24, a34 = 0.7, -0.3, 0.9, 0.4, -0.6, 0.2
    A = [0 a12 a13 a14; -a12 0 a23 a24; -a13 -a23 0 a34; -a14 -a24 -a34 0]
    @test wick_pfaffian(A) ≈ a12 * a34 - a13 * a24 + a14 * a23 atol = 1e-12
    # the defining identity Pf(A)² = det(A) (independent construction)
    @test wick_pfaffian(A)^2 ≈ det(A) atol = 1e-10
    # a larger random antisymmetric 6×6
    R = [i < j ? (0.13 * i - 0.07 * j + 0.5) : 0.0 for i in 1:6, j in 1:6]
    Aanti = R - R'
    @test wick_pfaffian(Aanti)^2 ≈ det(Aanti) atol = 1e-8
    # odd dimension ⇒ Pfaffian 0
    @test wick_pfaffian(zeros(3, 3)) == 0.0
end

@testset "bosonic Wick: the permanent, contrasted with the fermion determinant" begin
    using AbstractQAtlas: wick_permanent, wick_contraction
    n = 0.8
    G = fill(n, 1, 1)                       # single mode, ⟨a†a⟩ = n
    # bosons: ⟨a†a†aa⟩ = perm = 2n² = ⟨n(n−1)⟩ (super-Poissonian Gaussian boson)
    @test wick_permanent(G, [1, 1], [1, 1]) ≈ 2 * n^2 atol = 1e-12
    # fermions: ⟨c†c†cc⟩ = det = 0 (Pauli exclusion) — same indices, opposite statistics
    @test wick_contraction(G, [1, 1], [1, 1]) ≈ 0.0 atol = 1e-14
    # 2×2 permanent perm([[a,b],[c,d]]) = ad + bc
    M = [0.4 0.6; 0.3 0.9]
    @test wick_permanent(M, [1, 2], [1, 2]) ≈ 0.4 * 0.9 + 0.6 * 0.3 atol = 1e-12
    @test_throws ErrorException wick_permanent(M, [1, 2], [1])
end

@testset "Bloch–De Dominicis thermal contractions (Fermi / Bose occupations)" begin
    ε, β = 0.7, 1.5
    @test check(FermiDiracContraction(); n=1 / (exp(β * ε) + 1), ε=ε, β=β, atol=1e-12)
    @test check(BoseEinsteinContraction(); n=1 / (exp(β * ε) - 1), ε=ε, β=β, atol=1e-12)
    # β↔T convention, and the physical limits
    @test solve(FermiDiracContraction(), Val(:n); ε=ε, T=1 / β) ≈ 1 / (exp(β * ε) + 1)
    @test solve(FermiDiracContraction(), Val(:n); ε=0.0, β=β) ≈ 0.5      # ε=0 ⇒ half-filling
    @test solve(FermiDiracContraction(), Val(:n); ε=100.0, β=β) ≈ 0.0 atol = 1e-9  # empty high above E_F
end
