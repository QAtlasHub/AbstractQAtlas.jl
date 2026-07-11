# relations/wick.jl — Wick's theorem for Gaussian fermion states.
#
# Number-conserving case only in v0.1: every 2n-point function of a
# Gaussian (free-fermion, thermal or ground-state) density matrix reduces
# to a determinant of the 2-point matrix  G_ij = ⟨c†_i c_j⟩.
# The anomalous / BdG generalization (Pfaffian of the full covariance,
# ⟨cc⟩ ≠ 0) is a tracked follow-up.

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
