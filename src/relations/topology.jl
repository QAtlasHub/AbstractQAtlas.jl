# relations/topology.jl — standard topological invariants on Bloch maps.
#
# Generic numerical *definitions* operating on user-supplied Bloch
# functions — no model types, no atlas values.  Strictly textbook-level:
# 1D two-band winding, lattice (Fukui–Hatsugai–Suzuki) Chern number, and
# the TKNN quantization statement.
#
# References: Thouless–Kohmoto–Nightingale–den Nijs, PRL 49, 405 (1982);
# Fukui, Hatsugai & Suzuki, J. Phys. Soc. Jpn. 74, 1674 (2005).

"""
    winding_number(dvec::Function; nk::Int=1001) -> Int

Winding number of a planar map `k ∈ [0, 2π) ↦ dvec(k) = (d_x, d_y)`
around the origin — the 1D two-band invariant (e.g. SSH:
`d(k) = (v + w cos k, w sin k)` winds once for `v < w`, zero times for
`v > w`).

Computed by accumulating the exterior-angle increments
`Δθ = atan(d₁×d₂, d₁⋅d₂)` between successive grid points and rounding
the total to an integer (exact for a polygon avoiding the origin).

Throws when the map is not resolved: `|d(k)|` must stay larger than
twice the local polygon step everywhere, otherwise the curve passes
within a grid step of the origin — a gap closing (invariant undefined)
or an under-resolved map (increase `nk`).  A grid-point check alone is
NOT sufficient: a curve through the origin *between* samples still
yields an integer polygon winding, silently wrong.
"""
function winding_number(dvec::Function; nk::Int=1001)
    ks = range(0, 2π; length=nk + 1)  # closed loop: last point == first
    pts = [dvec(k) for k in ks]
    for m in 1:nk
        (x1, y1) = pts[m]
        (x2, y2) = pts[m + 1]
        r = hypot(x1, y1)
        step = hypot(x2 - x1, y2 - y1)
        r > 1e-12 || error("winding_number: |d(k)| ≈ 0 at k=$(ks[m]) — gapless")
        r > 2 * step || error(
            "winding_number: |d(k)| at k=$(ks[m]) is within two grid steps of " *
            "the origin — gap closing or under-resolved map; increase nk",
        )
    end
    total = 0.0
    for m in 1:nk
        (x1, y1) = pts[m]
        (x2, y2) = pts[m + 1]
        total += atan(x1 * y2 - y1 * x2, x1 * x2 + y1 * y2)
    end
    return round(Int, total / (2π))
end
export winding_number

"""
    chern_number(hk::Function, nbands::Int; nk::Int=24) -> Int

Chern number of the lowest `nbands` bands of a Bloch Hamiltonian
`hk(kx, ky) -> AbstractMatrix` (Hermitian), via the
Fukui–Hatsugai–Suzuki lattice field-strength method on an `nk × nk`
Brillouin-zone grid.

FHS is gauge-invariant by construction and returns the *exact* integer
already on coarse grids (provided the grid resolves the gap) — that is
the method's selling point, and why the result is `round`ed to `Int`
with a large-deviation guard rather than reported as a float.

Throws if the spectral gap between band `nbands` and `nbands + 1`
(numerically) closes anywhere on the grid.
"""
function chern_number(hk::Function, nbands::Int; nk::Int=24)
    # Occupied-frame eigenvectors on the grid.
    frames = Matrix{Matrix{ComplexF64}}(undef, nk, nk)
    ks = [2π * (m - 1) / nk for m in 1:nk]
    for (ix, kx) in enumerate(ks), (iy, ky) in enumerate(ks)
        H = Matrix{ComplexF64}(hk(kx, ky))
        vals, vecs = eigen(Hermitian(H))
        gap = vals[nbands + 1] - vals[nbands]
        gap > 1e-10 || error(
            "chern_number: gap between bands $nbands and $(nbands + 1) closes " *
            "at (kx, ky) = ($kx, $ky) (gap = $gap) — invariant undefined",
        )
        frames[ix, iy] = vecs[:, 1:nbands]
    end
    # Link variables U_μ(k) = det(V†(k) V(k+μ̂)) and plaquette field
    # strength F = arg(U₁ U₂ U₁⁻¹ U₂⁻¹), principal branch.
    wrap(i) = mod1(i, nk)
    total = 0.0
    for ix in 1:nk, iy in 1:nk
        V00 = frames[ix, iy]
        V10 = frames[wrap(ix + 1), iy]
        V11 = frames[wrap(ix + 1), wrap(iy + 1)]
        V01 = frames[ix, wrap(iy + 1)]
        u1 = det(V00' * V10)
        u2 = det(V10' * V11)
        u3 = det(V11' * V01)
        u4 = det(V01' * V00)
        total += angle(u1 * u2 * u3 * u4)
    end
    c = total / (2π)
    abs(c - round(c)) < 1e-6 || error(
        "chern_number: lattice sum $(c) is not integral — grid too coarse " *
        "for this gap; increase nk",
    )
    return round(Int, c)
end
export chern_number

"""
    TKNN <: AbstractRelation

The TKNN quantization statement: the zero-temperature Hall conductivity
of a gapped 2D band insulator is `σ_xy = C · e²/h`, with `C` the total
Chern number of the occupied bands.

Variables (in units of `e²/h`): `σxy`, `C`.
"""
@relation :topology TKNN(σxy, C) = σxy - C
