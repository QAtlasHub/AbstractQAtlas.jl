# relations/wick.jl — Wick's theorem and its finite-temperature
# (Bloch–De Dominicis) generalization for Gaussian states.
#
# Every 2n-point function of a Gaussian (free-fermion / free-boson)
# density matrix — ground state OR thermal — factorizes into a sum over
# all pairings of two-point contractions:
#   * number-conserving fermions → DETERMINANT of G_ij = ⟨c†_i c_j⟩
#     (Wick, Phys. Rev. 80, 268 (1950); at finite T this is the
#     Bloch–De Dominicis theorem, Nucl. Phys. 7, 459 (1958), with the
#     THERMAL contraction G_ij = δ_ij/(e^{βε_i}+1));
#   * general (BdG / paired, ⟨cc⟩ ≠ 0) fermions → PFAFFIAN of the
#     antisymmetric contraction matrix;
#   * bosons → PERMANENT of ⟨a†_i a_j⟩ (same pairings, all + signs; the
#     thermal contraction is the Bose factor 1/(e^{βε}−1)).
# The state never enters — only its two-point data does.  That IS the
# theorem: Gaussian states are fixed by their contractions.

"""
    wick_contraction(G::AbstractMatrix, cr::AbstractVector{<:Integer},
                     an::AbstractVector{<:Integer}) -> Number

Wick's theorem for a number-conserving Gaussian fermion state with
2-point matrix `G[i, j] = ⟨c†_i c_j⟩`:

`⟨c†_{cr[1]} ⋯ c†_{cr[n]} c_{an[n]} ⋯ c_{an[1]}⟩ = det(M)`,
`M[p, q] = G[cr[p], an[q]]`.

**Ordering convention**: creation operators appear left-to-right as
`cr[1], …, cr[n]`; annihilation operators appear in *reversed* order
`an[n], …, an[1]` (the nested ordering, so `cr = [i], an = [j]` gives
`⟨c†_i c_j⟩ = G[i, j]` and `cr = an = [i, j]` gives
`⟨c†_i c†_j c_j c_i⟩ = ⟨n_i n_j⟩` for `i ≠ j`).

The state itself never enters — only `G` does.  That *is* Wick's
theorem: Gaussian states are fully determined by their 2-point data.
"""
function wick_contraction(
    G::AbstractMatrix, cr::AbstractVector{<:Integer}, an::AbstractVector{<:Integer}
)
    length(cr) == length(an) ||
        error("particle-number conservation requires equal counts of c† and c")
    n = length(cr)
    M = [G[cr[p], an[q]] for p in 1:n, q in 1:n]
    return det(M)
end
export wick_contraction

"""
    wick_density_correlation(G::AbstractMatrix, i::Integer, j::Integer) -> Number

Density–density correlation of a number-conserving Gaussian fermion
state from its 2-point matrix:

`⟨n_i n_j⟩ = G_ii G_jj + G_ij (δ_ij − G_ji)`.

Derivation: for `i ≠ j`, Wick gives
`⟨c†_i c†_j c_j c_i⟩ = G_ii G_jj − G_ij G_ji`; for `i = j`, `n_i² = n_i`
gives `G_ii`.  Both cases are captured by the single formula above.
"""
function wick_density_correlation(G::AbstractMatrix, i::Integer, j::Integer)
    return G[i, i] * G[j, j] + G[i, j] * ((i == j ? 1 : 0) - G[j, i])
end
export wick_density_correlation

# ─── Anomalous / BdG Wick: the Pfaffian (issue #5) ──────────────────────

"""
    wick_pfaffian(A::AbstractMatrix) -> Number

The Pfaffian of a `2n × 2n` antisymmetric matrix `A` — the
number-**non**-conserving (BdG / paired) form of Wick's theorem: for a
general Gaussian state the `2n`-point Majorana correlation is the sum over
all pairings, `Pf(A)`, where `A_ab = ⟨γ_a γ_b⟩` is the antisymmetric
contraction matrix.  Reduces to the [`wick_contraction`](@ref) determinant
when the state is number-conserving (no anomalous `⟨cc⟩` pairing), since
`Pf(A)² = det(A)`.

Computed by the exact cofactor recursion `Pf(A) = Σ_{j≥2} (−1)^j A_{1j}
Pf(A_{1̂ĵ})` (small `n`; `Pf = 0` for odd dimension, `1` for the empty
matrix).
"""
function wick_pfaffian(A::AbstractMatrix)
    n = size(A, 1)
    n == size(A, 2) || error("wick_pfaffian: matrix must be square")
    isodd(n) && return zero(eltype(A))
    n == 0 && return one(eltype(A))
    pf = zero(eltype(A))
    for j in 2:n
        rest = [k for k in 2:n if k != j]
        pf += (-1)^j * A[1, j] * wick_pfaffian(A[rest, rest])
    end
    return pf
end
export wick_pfaffian

# ─── Bosonic Wick: the permanent (issue #9) ─────────────────────────────

"""
    wick_permanent(G::AbstractMatrix, cr::AbstractVector{<:Integer},
                   an::AbstractVector{<:Integer}) -> Number

Wick's theorem for a Gaussian **boson** state with 2-point matrix
`G[i, j] = ⟨a†_i a_j⟩`: the same sum over pairings as the fermionic
[`wick_contraction`](@ref) but with **all `+` signs**, i.e. the PERMANENT
instead of the determinant,

`⟨a†_{cr[1]} ⋯ a†_{cr[n]} a_{an[n]} ⋯ a_{an[1]}⟩ = perm(M)`,
`M[p, q] = G[cr[p], an[q]]`.

(The Bose symmetry replaces the fermionic antisymmetric sign; at finite T
the contraction is the Bose factor.)
"""
function wick_permanent(
    G::AbstractMatrix, cr::AbstractVector{<:Integer}, an::AbstractVector{<:Integer}
)
    length(cr) == length(an) ||
        error("particle-number conservation requires equal counts of a† and a")
    n = length(cr)
    M = [G[cr[p], an[q]] for p in 1:n, q in 1:n]
    return _permanent(M)
end
export wick_permanent

# permanent by cofactor expansion along the first row (small n)
function _permanent(M::AbstractMatrix)
    n = size(M, 1)
    n == 0 && return one(eltype(M))
    n == 1 && return M[1, 1]
    total = zero(eltype(M))
    for j in 1:n
        cols = [c for c in 1:n if c != j]
        total += M[1, j] * _permanent(M[2:n, cols])
    end
    return total
end

# ─── Bloch–De Dominicis thermal contractions ────────────────────────────

"""
    FermiDiracContraction <: AbstractRelation

The finite-temperature two-point contraction of the Bloch–De Dominicis
theorem (Bloch & De Dominicis, Nucl. Phys. 7, 459 (1958)) for fermions —
the mode occupation that seeds the thermal Wick determinant/Pfaffian,

`⟨c†_ε c_ε⟩ = n_F(ε) = 1/(e^{βε} + 1)`.

Variables: `n`, `ε`, `β` (or `T`).
"""
@relation :wick FermiDiracContraction(n, ε, β) = n - 1 / (exp(β * ε) + 1)

"""
    BoseEinsteinContraction <: AbstractRelation

The finite-temperature two-point contraction for bosons — the mode
occupation seeding the thermal Wick **permanent**,

`⟨a†_ε a_ε⟩ = n_B(ε) = 1/(e^{βε} − 1)`   (`ε > 0`).

Variables: `n`, `ε`, `β` (or `T`).
"""
@relation :wick BoseEinsteinContraction(n, ε, β) = n - 1 / (exp(β * ε) - 1)
