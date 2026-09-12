# relations/quantity_links.jl — the relation → quantity map.
#
# Declares, for each relation whose subject is a named quantity, which
# `AbstractQuantity` TYPES it directly constrains (parametric families use
# the UnionAll, so `Susceptibility` matches any `Susceptibility{I}`).  This
# is the machine link that turns the registry into a queryable web:
# `relations_constraining(q)` reads it in reverse.  Included after the
# physics submodules are re-exported, so every relation + quantity name
# resolves at the top level.  Relations that constrain parameters or
# exponents rather than named quantities (scaling laws, Maxwell relations,
# the uncertainty relation) keep the default `quantities(rel) == ()`.

# `also_constrains` hints — for TYPE-KEYED relations, a quantity that a relation
# constrains through a SUPPLIED variance/derivative (not a typed subject) is invisible
# to auto-derivation; declared here so `quantities` (= typed subjects ∪ these) keeps
# the physics graph complete (e.g. Var(M) → Magnetization in the susceptibility FDT):
also_constrains(::SusceptibilityFDT) = (Magnetization,)
also_constrains(::SpecificHeatFDT) = (Energy,)
also_constrains(::SpecificHeatFromEntropy) = (ThermalEntropy,)
also_constrains(::MottFormula) = (Conductivity,)   # transport: dlnσ/dε → Conductivity
also_constrains(::MicrocanonicalTemperature) = (ThermalEntropy, Energy)  # β = ∂S/∂E
also_constrains(::EntropyResponse) = (FreeEnergy,)     # fundamental: S = −∂F/∂T
also_constrains(::GibbsHelmholtz) = (FreeEnergy,)      # fundamental: U = ∂(βF)/∂β
also_constrains(::MagnetizationResponse) = (FreeEnergy,)     # fundamental: M = −∂F/∂h
also_constrains(::SusceptibilityResponse) = (Magnetization,)  # fundamental: χ = ∂M/∂h
also_constrains(::ParticleNumberResponse) = (GrandPotential,)  # grand-canonical: N = −∂Ω/∂μ
also_constrains(::StaticFromDynamicalStructureFactor) = (DynamicalStructureFactor,)  # Sq = ∫S(q,ω)dω/2π (supplied)
also_constrains(::ChernFromBerryCurvature) = (BerryCurvature,)  # topology: C = ∫Ω d²k/2π (supplied integral)
also_constrains(::CFTEntanglementSlope) = (VonNeumannEntropy,)  # entanglement: dS/d(ln ℓ) (supplied derivative)
also_constrains(::DynamicalScaling) = (MassGap,)            # scaling: d(lnΔ)/d(lnξ) (supplied)
# The infinite-randomness block keys on the REDUCTIONS, not the clean quantity.
# `Typical` and `DisorderAveraged` are separate quantities here precisely so a bag
# cannot mix them (core/quantities.jl), and every relation below is a statement
# about one of them — a bare key would put the clean value in the same slot.
# `TypicalCorrelationLength` is the sharp case: it exists to say the two lengths
# differ, so keying it on `CorrelationLength` names the thing it denies.
also_constrains(::ActivatedDynamicalScaling) = (MassGap, Typical{MassGap})  # ψ is the TYPICAL gap's
also_constrains(::GriffithsExponentDivergence) = (DynamicalExponent,)  # an ensemble exponent, not a reduction
also_constrains(::GriffithsSusceptibility) = (DisorderAveraged{Susceptibility},)  # rare regions ⇒ the average
also_constrains(::GriffithsSpecificHeat) = (DisorderAveraged{SpecificHeat},)      # same
function also_constrains(::TypicalCorrelationLength)
    return (Typical{CorrelationLength}, DisorderAveraged{CorrelationLength})
end
# Appendix A: the size-space and dynamics forms carry the same hidden subjects
# as their ξ-space partners above, with the same reduction split.
also_constrains(::ActivatedFiniteSizeScaling) = (MassGap, Typical{MassGap})
also_constrains(::ConventionalFiniteSizeEnergy) = (MassGap,)
also_constrains(::OrderParameterDimension) = (SpontaneousMagnetization,)
also_constrains(::CriticalAutocorrelation) = (DynamicalCorrelation,)
also_constrains(::CriticalQuantumSusceptibility) = (Susceptibility,)
also_constrains(::CriticalQuantumSpecificHeat) = (SpecificHeat,)
also_constrains(::ActivatedAutocorrelation) = (DisorderAveraged{DynamicalCorrelation},)
also_constrains(::ActivatedSusceptibility) = (DisorderAveraged{Susceptibility},)
also_constrains(::ActivatedSpecificHeat) = (DisorderAveraged{SpecificHeat},)
also_constrains(::GriffithsAutocorrelation) = (DisorderAveraged{DynamicalCorrelation},)
also_constrains(::ConventionalFieldSusceptibility) = (Susceptibility,)
also_constrains(::ConventionalFieldSpecificHeat) = (SpecificHeat,)
also_constrains(::OrderedGriffithsEnergyScale) = (MassGap,)
# HeatCapacityDifference GAINS ThermalExpansionCoefficient + IsothermalCompressibility
# (α, κT are now typed subjects → auto).  LinearResponseFDT now types `β::InverseTemperature`
# (bag-visible, like Jarzynski/Crooks) but has no quantity subject; the 4 Maxwell relations
# stay fully symbol-keyed (generic / pure-derivative — no named quantity subject).

# ── Correlations / Green's functions ──
# NB: Dyson, SpectralFromGreens and the whole Keldysh block below are now
# TYPE-KEYED (`Name(x::Quantity, …)` in spectral.jl / keldysh.jl), so their
# `quantities` is auto-derived from the declaration — no hand-link here.  The
# remaining entries are legacy symbol-keyed relations awaiting migration.
quantities(::SpectralSumRule) = (SpectralFunction,)
quantities(::FSumRule) = (DynamicalStructureFactor,)   # first moment of S(q,ω) (supplied)
quantities(::DetailedBalance) = (DynamicalStructureFactor,)
quantities(::DynamicalFDT) = (DynamicalStructureFactor, DynamicalSusceptibility)
quantities(::ResponseRealityReal) = (DynamicalSusceptibility,)   # reality: Re χ even under ω→−ω
quantities(::ResponseRealityImag) = (DynamicalSusceptibility,)   # reality: Im χ odd under ω→−ω
# CorrelationLengthGap (ξ::CorrelationLength, v::Velocity, Δ::MassGap), NMRExponent
# (θ_NMR::NMRRelaxationExponent, Δ_op::ScalingDimension), FiniteSizeGap, and
# StaticFromDynamicalStructureFactor (Sq::StaticStructureFactor, + also_constrains
# DynamicalStructureFactor for the supplied sqw_integral) are now type-keyed — `quantities`
# is auto-derived, no hand-link.  (The old NMRExponent hand-link named NMRSpinRelaxationRate,
# which is not even a variable of the relation.)

# ── Keldysh RAK structure & fluctuation–dissipation: now type-keyed (keldysh.jl),
#    `quantities` auto-derived — see the note above. ──

# ── Transport: now type-keyed (transport.jl) with CONCRETE component types
#    (Conductivity{(:x,:x)} etc.), `quantities` auto-derived. OnsagerReciprocity
#    (generic L_μν) and IoffeRegel (dimensionless k_Fℓ) stay symbol-keyed — no
#    named quantity subject, so `quantities` == () by default.
#    MottFormula's dlnσ/dε → Conductivity association is kept via `also_constrains`
#    above (it is a supplied derivative, not a typed subject); RighiLeduc also
#    constrains Conductivity now that σxy is typed. ──

# ── Quantum information & entanglement ──
quantities(::RenyiTwoPurity) = (RenyiEntropy, Purity)
quantities(::RenyiEntropyMoment) = (RenyiEntropy,)
quantities(::TsallisEntropyMoment) = (TsallisEntropy,)
quantities(::MutualInformationDefinition) = (MutualInformation,)
quantities(::ConditionalEntropyDefinition) = (ConditionalEntropy,)
quantities(::MeasurementEntropyIncrease) = (MeasurementEntropy,)
quantities(::MeasurementEntropyRelative) = (MeasurementEntropy, RelativeEntropy)
quantities(::MarkovEntropyDefinition) = (MarkovEntropy,)
quantities(::ConcurrenceTangle) = (Concurrence, Tangle)
quantities(::Monogamy) = (Tangle,)
quantities(::ThreeTangleDefinition) = (ThreeTangle, Tangle)
quantities(::TripartiteInformationDefinition) = (TripartiteInformation,)
quantities(::KitaevPreskillTEE) = (TopologicalEntanglementEntropy,)
quantities(::EntropyNonNegativity) = (VonNeumannEntropy,)
quantities(::MaxEntropyBound) = (VonNeumannEntropy,)
quantities(::Subadditivity) = (VonNeumannEntropy,)
quantities(::ArakiLieb) = (VonNeumannEntropy,)
quantities(::StrongSubadditivity) = (VonNeumannEntropy,)
quantities(::WeakMonotonicity) = (VonNeumannEntropy,)
quantities(::EntropyMixingConcavity) = (VonNeumannEntropy,)
quantities(::HolevoMixingBound) = (VonNeumannEntropy,)
quantities(::RenyiMonotonicity) = (RenyiEntropy,)
quantities(::RelativeEntropyNonNegativity) = (RelativeEntropy,)

# ── Disorder statistics ──
# Both are fully symbol-keyed: `X_typ`/`X_avg` and `F_quenched`/`F_annealed` are
# REDUCTIONS over an ensemble, and no single quantity names a reduction, so there
# is no slot to type.  Linked by hand for the same reason the entropy
# inequalities above are — otherwise `relations_constraining(FreeEnergy)` cannot
# find the annealed bound, and `Typical`/`DisorderAveraged` are orphans of the
# one relation that is about them.
quantities(::TypicalBelowAverage) = (Typical, DisorderAveraged)
quantities(::AnnealedFreeEnergyBound) = (FreeEnergy,)

# ── Quantum-mechanical foundations ──
# VirialTheorem is type-keyed (quantum.jl), `quantities` auto-derived. The Ehrenfest /
# Hellmann–Feynman / uncertainty / Lieb–Robinson relations stay symbol-keyed (generic
# operators / derivatives). EnergyVarianceEigenstate stays symbol-keyed — its subject
# EnergyVariance enters only through the moment combination ⟨H²⟩ − E², not a slot.
quantities(::EnergyVarianceEigenstate) = (EnergyVariance,)

# ── Topology ──
# TKNN + ChernFromBerryCurvature are type-keyed (topology.jl), `quantities`
# auto-derived (TKNN → Conductivity, ChernNumber). BulkBoundary stays symbol-keyed
# (generic bulk invariant ν — Chern / winding / ℤ₂, no single named type).
