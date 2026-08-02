# relations/fluctuation.jl — nonequilibrium work fluctuation theorems.
#
# The universal work ↔ free-energy relations a driven (nonequilibrium) calculation
# checks against: Jarzynski's equality tying the exponential work average to the
# EQUILIBRIUM free-energy difference for ANY protocol, its second-law corollary
# ⟨W⟩ ≥ ΔF, and the Crooks theorem relating the forward/reverse work distributions.
# Domain tag :fluctuation.  The averages ⟨e^{−βW}⟩ / ⟨W⟩ and the distribution ratio
# are caller-supplied aggregates over the work statistics (the sampling is a functional
# step — cf. the supplied-average conventions in ensembles.jl).
#
# References (doiget-verified, docs/references.bib): Jarzynski, Phys. Rev. Lett. 78,
# 2690 (1997); Crooks, [Crooks1999](@cite).

"""
    JarzynskiEquality <: AbstractRelation

Jarzynski's nonequilibrium equality (Jarzynski, [Jarzynski1997](@cite)):
the exponential average of the work `W` over many realizations of ANY protocol driving
the system between two equilibrium states equals the exponentiated equilibrium
free-energy difference, however far from equilibrium the driving is,

`⟨e^{−βW}⟩ = e^{−β ΔF}`,   `ΔF = F_B − F_A`.

Supplied-average convention: `exp_work = ⟨e^{−βW}⟩`.

Variables: `exp_work` = `⟨e^{−βW}⟩`, `ΔF`, `β` (or `T`).
"""
@relation :fluctuation JarzynskiEquality(exp_work, ΔF, β::InverseTemperature) =
    exp_work - exp(-β * ΔF)

"""
    JarzynskiSecondLaw <: AbstractInequality

The second law as the Jensen-inequality corollary of [`JarzynskiEquality`](@ref)
(`⟨e^{−βW}⟩ ≥ e^{−β⟨W⟩}` by convexity ⇒ `e^{−βΔF} ≥ e^{−β⟨W⟩}`): the average work done
on the system cannot be less than the free-energy difference,

`⟨W⟩ ≥ ΔF`

(slack `W_avg − ΔF` = the dissipated work `W_diss ≥ 0`).  Saturated by a quasistatic
(reversible) protocol; a strictly positive slack measures irreversibility.

Variables: `ΔF` (the bounded one), `W_avg` = `⟨W⟩`.
"""
@bound :fluctuation JarzynskiSecondLaw(ΔF <= W_avg)

"""
    CrooksFluctuationTheorem <: AbstractRelation

The Crooks fluctuation theorem (Crooks, [Crooks1999](@cite)): the forward and
time-reversed work distributions of a driven process obey

`P_F(W) / P_R(−W) = e^{β (W − ΔF)}`,

crossing at `W = ΔF` (`ratio = 1`) and integrating over `W` to [`JarzynskiEquality`](@ref).
Supplied-ratio convention: `ratio = P_F(W) / P_R(−W)`.

Variables: `ratio` = `P_F(W)/P_R(−W)`, `W`, `ΔF`, `β` (or `T`).
"""
@relation :fluctuation CrooksFluctuationTheorem(ratio, W, ΔF, β::InverseTemperature) =
    ratio - exp(β * (W - ΔF))
