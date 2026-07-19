# The functional-evaluation interface (scope-line #14 seam): the generic verbs are owned
# here (like `fetch`), the numerics belong to the functional sibling.  Two things to
# check: the "not implemented" fallback errors informatively, and a MINIMAL sibling
# implementation wires the verbs end-to-end to the supplied-integral relations.

using AbstractQAtlas
using AbstractQAtlas: principal_value_hilbert, spectral_moment
using Test

struct _NoImplResponse <: AbstractResponse end

@testset "evaluation seam: the functional-sibling fallback errors informatively" begin
    r = _NoImplResponse()
    @test_throws ErrorException principal_value_hilbert(r, 0.5)
    @test_throws ErrorException spectral_moment(r, 0)
    e = try
        spectral_moment(r, 0)
    catch e
        e
    end
    @test occursin("functional-sibling", e.msg) && occursin("_NoImplResponse", e.msg)
end

# a MINIMAL functional-sibling implementation — a dense (ω, values) grid + trapezoid
# quadrature.  The real numerics live in the functional package; this is just enough to
# prove the interface wires end-to-end to the supplied-integral relations.
struct _GridResponse <: AbstractResponse
    ω::Vector{Float64}
    f::Vector{Float64}
end
function AbstractQAtlas.spectral_moment(r::_GridResponse, n::Integer)
    g = (r.ω .^ n) .* r.f
    return sum((g[1:(end - 1)] .+ g[2:end]) ./ 2 .* diff(r.ω))       # ∫ ωⁿ f dω (trapezoid)
end

@testset "seam → relation: a grid response gates the sum rules turnkey" begin
    ω = collect(range(-150.0, 150.0; length=300_001))
    η = 0.05
    # 0th moment of a Lorentzian A(ω) feeds the spectral sum rule ∫A = 1
    A = (1 / π) .* η ./ ((ω .- 1.3) .^ 2 .+ η^2)
    Arep = _GridResponse(ω, A)
    @test isapprox(spectral_moment(Arep, 0), 1.0; atol=1e-3)
    @test check(SpectralSumRule(); spectral_integral=spectral_moment(Arep, 0), atol=1e-3)
    # 1st moment of a single-mode S(q,ω) (Gaussian at ω0 = q²/2m) feeds the f-sum ∫ωS = q²/2m
    q, m = 1.4, 1.0
    ω0 = q^2 / (2m)
    S = (1 / (η * sqrt(2π))) .* exp.(-((ω .- ω0) .^ 2) ./ (2η^2))
    Srep = _GridResponse(ω, S)
    @test isapprox(spectral_moment(Srep, 1), ω0; atol=1e-3)
    @test check(FSumRule(); first_moment=spectral_moment(Srep, 1), q=q, m=m, atol=1e-3)
end
