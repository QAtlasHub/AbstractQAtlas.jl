# Reading a central charge off measured entanglement, and what a convention
# declaration buys.
#
# The numbers are not pasted in from anywhere: a critical free-fermion chain's
# block correlation matrix is a closed form, so the entropies below are computed
# here and the central charge that comes out is a measurement, not a fixture. Its
# distance from the exact `c = 1` is the finite-block correction, which is why
# `test/core/test_examples.jl` checks it to 1e-3 and not to machine precision.
#
# Wrapped in a module because `test/core/test_examples.jl` includes this file, and
# every test file is included into one namespace: `block_entropy` and a couple of
# block sizes at top level would be in every other test file's way.
#
# Run: julia --project=. examples/critical_entanglement.jl

module CriticalEntanglementExample

using AbstractQAtlas
using LinearAlgebra: eigvals, Symmetric

export block_entropy, central_charges

"""
    correlation_matrix(ℓ) -> Symmetric

The `ℓ x ℓ` block of the half-filled critical free-fermion chain's two-point
function, `C_jk = sin(π(j-k)/2) / (π(j-k))`, `C_jj = 1/2`.
"""
function correlation_matrix(ℓ::Integer)
    return Symmetric([
        j == k ? 0.5 : sin(π * (j - k) / 2) / (π * (j - k)) for j in 1:ℓ, k in 1:ℓ
    ])
end

"Entanglement entropy of a block of `ℓ` sites, in nats."
function block_entropy(ℓ::Integer)
    return free_fermion_entanglement_entropy(eigvals(correlation_matrix(ℓ)))
end

"""
    central_charges(small = 16, large = 64) -> NamedTuple

`c` read off the same measurement entered three ways: in nats, in bits with
nothing declared, and in bits declared.

Two block sizes, because the charge is not falsifiable from one: the logarithmic
forms each carry a non-universal constant that can absorb any `c`, and the slope
through two blocks carries none. The middle entry is the failure the convention
layer exists for, and it is wrong by exactly the base nobody mentioned.
"""
function central_charges(small::Integer=16, large::Integer=64)
    a, b = Region((1:small)...), Region((1:large)...)
    span = log(large) - log(small)
    slope(bg) = (bg[entanglement_entropy(b)] - bg[entanglement_entropy(a)]) / span
    pack(f) = (entanglement_entropy(a) => f(small), entanglement_entropy(b) => f(large))
    # Read out of the bag, so a declared convention reaches the slope through it.
    c(bg) = derive_crosschecked(:c; dS_dlogℓ=slope(bg), ncuts=2)
    in_bits(ℓ) = block_entropy(ℓ) / log(2)
    return (
        nats=c(bag(pack(block_entropy)...)),
        bits_undeclared=c(bag(pack(in_bits)...)),
        bits_declared=c(
            bag(conventions(AbstractEntanglementMeasure => Bits), pack(in_bits)...)
        ),
    )
end

end # module CriticalEntanglementExample

if abspath(PROGRAM_FILE) == @__FILE__
    let cs = CriticalEntanglementExample.central_charges()
        println("c from blocks of 16 and 64 sites, exact value 1:")
        println("  nats, as the calculation produced them : ", round(cs.nats; digits=4))
        println(
            "  the same numbers in bits, undeclared   : ",
            round(cs.bits_undeclared; digits=4),
        )
        println(
            "      which is 1/ln2 = ",
            round(1 / log(2); digits=4),
            ", the base it was never told",
        )
        println(
            "  the same numbers in bits, declared     : ", round(cs.bits_declared; digits=4)
        )
    end
end
