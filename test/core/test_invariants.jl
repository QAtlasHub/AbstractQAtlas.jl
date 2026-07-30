# Structural invariants that must hold across the WHOLE quantity vocabulary
# — enumerated by reflection, so a newly-added quantity is checked
# automatically (and a one-sided declaration, as `fourier_pair` once had,
# fails CI).  These are the universal laws the tensor traits promise, not
# per-quantity facts.

using AbstractQAtlas
using AbstractQAtlas:
    tensor_rank,
    index_spaces,
    indices,
    frequency_arguments,
    response_order,
    fourier_conjugate_quantity,
    fourier_pair
using InteractiveUtils: subtypes

# every concrete leaf of the AbstractQuantity tree
function _leaves(T, acc=Type[])
    for S in subtypes(T)
        isabstracttype(S) ? _leaves(S, acc) : push!(acc, S)
    end
    return acc
end

# a representative instance: try the common constructor shapes, plus the
# `ThermalAverage` wrapper (which delegates its tensor traits to the
# wrapped quantity).
function _representative(T)
    T === ThermalAverage && return ThermalAverage(Susceptibility(:x, :y), Canonical(1.0))
    # `(2,)` last: a numeric shape would otherwise shadow a symbol-taking constructor
    # that happens to accept an Int, and the first success wins. It is here for the
    # quantities whose ORDER is data rather than a type parameter — `RenyiEntropy(α)`
    # is the one, and α = 1 is refused by construction (it is the von Neumann limit),
    # so the representative must not be 1.
    for a in ((), (:x, :y), (:x, :y, :z), (:x,), (:natural,), (2,))
        try
            return T(a...)
        catch
        end
    end
    return nothing
end

const _QUANTITY_LEAVES = _leaves(AbstractQAtlas.AbstractQuantity)

@testset "every concrete quantity is instantiable (reflection coverage)" begin
    @test !isempty(_QUANTITY_LEAVES)
    for T in _QUANTITY_LEAVES
        # a leaf the representative chain cannot build is a coverage gap —
        # add a constructor shape rather than silently skipping it.
        @test _representative(T) !== nothing
    end
end

@testset "tensor-trait invariants (length(index_spaces) == tensor_rank == length(indices))" begin
    for T in _QUANTITY_LEAVES
        q = _representative(T)
        q === nothing && continue
        r = tensor_rank(q)
        @test r >= 0
        @test length(index_spaces(q)) == r          # one index space per tensor slot
        idx = indices(q)
        isempty(idx) || @test length(idx) == r       # indices, when given, fill the rank
        @test frequency_arguments(q) >= 0
        @test response_order(q) >= 0
    end
end

@testset "response quantities: response_order == tensor_rank − 1" begin
    # a response tensor carries one output index + n field indices
    for q in (
        Susceptibility(:x, :y),
        Susceptibility(:x, :y, :z),
        DynamicalSusceptibility(:x, :y, :z),
        Conductivity(:x, :y),
        DynamicalConductivity(:x, :y, :z),
    )
        @test response_order(q) == tensor_rank(q) - 1
    end
end

@testset "Fourier conjugacy is an involution + fourier_pair is symmetric (all declared)" begin
    # every quantity that declares a conjugate must round-trip to its own
    # family, and pair symmetrically — generalising the per-pair checks so a
    # future one-sided fourier_conjugate_quantity declaration fails here.
    for T in _QUANTITY_LEAVES
        c = fourier_conjugate_quantity(T)
        c === nothing && continue
        back = fourier_conjugate_quantity(c)
        @test back !== nothing                       # the reverse must be declared
        @test T <: back                              # …and it round-trips to T's family
        a = _representative(T)
        a === nothing && continue
        # a matching conjugate instance (same indices if the family is parametric,
        # else the singleton)
        bi = try
            c(indices(a)...)
        catch
            try
                c()
            catch
                nothing
            end
        end
        bi === nothing && continue
        @test fourier_pair(a, bi) == fourier_pair(bi, a)   # symmetric
        @test fourier_pair(a, bi)                          # …and actually conjugate
    end
end
