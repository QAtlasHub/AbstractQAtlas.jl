# The dynamical-quantity inter-relationship graph: the operation edges,
# the chain from any derived quantity back to its source, and the tie
# from each edge to its exact pointwise relation (where one exists).

using AbstractQAtlas
using AbstractQAtlas: spectral_origin, spectral_chain, origin_relation

@testset "graph edges (which operation makes which quantity)" begin
    @test spectral_origin(RetardedGreensFunction()) == SpectralOrigin(SelfEnergy, :dyson)
    @test spectral_origin(SpectralFunction()) ==
        SpectralOrigin(RetardedGreensFunction, :neg_im_over_pi)
    @test spectral_origin(DensityOfStates()) ==
        SpectralOrigin(SpectralFunction, :bz_average)
    @test spectral_origin(DynamicalStructureFactor()) ==
        SpectralOrigin(DynamicalCorrelation, :spacetime_fourier)
    @test spectral_origin(NMRSpinRelaxationRate()) ==
        SpectralOrigin(DynamicalSusceptibility, :low_frequency_limit)
    # Kubo edge: an n-th order response is an n-time correlation — the
    # dynamical susceptibility routes to the SAME-ORDER correlation, so the
    # index tuple (hence the number of times) is preserved across the edge.
    @test spectral_origin(DynamicalSusceptibility(:x, :y)) ==
        SpectralOrigin(DynamicalCorrelation{(:x, :y)}, :kubo)          # linear ⟵ 2-point
    @test spectral_origin(DynamicalSusceptibility(:x, :y, :z)) ==
        SpectralOrigin(DynamicalCorrelation{(:x, :y, :z)}, :kubo)      # χ⁽²⁾(ω₁,ω₂) ⟵ 3-point
    # sources / off-graph quantities have no edge
    @test spectral_origin(SelfEnergy()) === nothing
    @test spectral_origin(DynamicalCorrelation(:x, :y)) === nothing
    @test spectral_origin(FreeEnergy()) === nothing
end

@testset "n-th order response ⟺ n-time correlation (Kubo edge is order-faithful)" begin
    # the physical invariant the user asked for: the Kubo correlation partner
    # of an order-n response has exactly n time differences = n frequencies,
    # matching the response's own frequency count — in BOTH channels (the
    # spin-response susceptibility and the current-response conductivity).
    using AbstractQAtlas: frequency_arguments, response_order
    for I in ((:x, :y), (:x, :y, :z), (:x, :x, :x, :x))
        χ = DynamicalSusceptibility(I...)
        corr = spectral_origin(χ).from            # the DynamicalCorrelation{I}
        n = response_order(χ)
        @test frequency_arguments(χ) == n                       # n-th order → n frequencies
        @test frequency_arguments(corr) == n                    # …and its correlation is n-time
        @test frequency_arguments(corr) == frequency_arguments(χ)

        σ = DynamicalConductivity(I...)
        jj = spectral_origin(σ).from              # the CurrentCorrelation{I}
        @test spectral_origin(σ) == SpectralOrigin(CurrentCorrelation{I}, :kubo)
        @test frequency_arguments(σ) == n
        @test frequency_arguments(jj) == frequency_arguments(σ)  # current corr is n-time too
        @test spectral_chain(σ) == [DynamicalConductivity{I}, CurrentCorrelation{I}]
    end
end

@testset "chains trace back to the source quantity" begin
    # the density of states is built from the self-energy: ρ ⟵ A ⟵ G^R ⟵ Σ
    @test spectral_chain(DensityOfStates()) ==
        [DensityOfStates, SpectralFunction, RetardedGreensFunction, SelfEnergy]
    @test spectral_chain(SpectralFunction()) ==
        [SpectralFunction, RetardedGreensFunction, SelfEnergy]
    @test spectral_chain(DynamicalStructureFactor()) ==
        [DynamicalStructureFactor, DynamicalCorrelation]
    # NMR ⟵ χ'' ⟵ (Kubo) correlation: 1/T₁ built from the correlation function
    @test spectral_chain(NMRSpinRelaxationRate()) ==
        [NMRSpinRelaxationRate, DynamicalSusceptibility, DynamicalCorrelation]
    @test spectral_chain(SelfEnergy()) == [SelfEnergy]           # source singleton
    # current channel mirrors the spin channel's fluctuation branch:
    # S^j(q,ω) ⟵ ⟨jj⟩(r,t) via the space-time Fourier transform
    @test spectral_origin(CurrentNoise(:x, :y)) ==
        SpectralOrigin(CurrentCorrelation, :spacetime_fourier)
    @test spectral_chain(CurrentNoise(:x, :y)) ==
        [CurrentNoise{(:x, :y)}, CurrentCorrelation]
end

@testset "edges tie to their exact pointwise relation (or nothing)" begin
    # single-(q,ω)-point operations point at the @relation that realizes them
    @test origin_relation(:dyson) isa Dyson
    @test origin_relation(:neg_im_over_pi) isa SpectralFromGreens
    # transform / BZ-sum / limit operations have no pointwise relation here
    # (their evaluation is the functional sibling's job — issue #14)
    @test origin_relation(:bz_average) === nothing
    @test origin_relation(:spacetime_fourier) === nothing
    @test origin_relation(:low_frequency_limit) === nothing
    @test origin_relation(:kubo) === nothing            # transform of a multi-time correlation
    # every graph edge whose via is pointwise resolves to a registered relation
    for q in (RetardedGreensFunction(), SpectralFunction())
        rel = origin_relation(spectral_origin(q).via)
        @test rel !== nothing
        @test rel in all_relations(; domain=:spectral)
    end
end

@testset "operation_scope: the definitional/functional line (issue #14)" begin
    using AbstractQAtlas: operation_scope
    # pointwise identities are definitional (owned here)
    @test operation_scope(:dyson) == :definitional
    @test operation_scope(:neg_im_over_pi) == :definitional
    # transforms / sums / limits are functional (the sibling evaluates them)
    for via in (:bz_average, :spacetime_fourier, :low_frequency_limit, :kubo)
        @test operation_scope(via) == :functional
    end
    # INVARIANT: the scope line is exactly origin_relation's split
    for via in (
        :dyson,
        :neg_im_over_pi,
        :bz_average,
        :spacetime_fourier,
        :low_frequency_limit,
        :kubo,
    )
        @test (operation_scope(via) == :definitional) == (origin_relation(via) !== nothing)
    end
    # every functional edge in the graph has no pointwise relation, and vice versa
    for q in (
        RetardedGreensFunction(),
        SpectralFunction(),
        DensityOfStates(),
        DynamicalStructureFactor(),
        DynamicalSusceptibility(:x, :y),
    )
        via = spectral_origin(q).via
        @test (operation_scope(via) == :functional) == (origin_relation(via) === nothing)
    end
end
