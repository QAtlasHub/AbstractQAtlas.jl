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

# ── Statistical mechanics / fundamental / ensemble (symbol-keyed, awaiting
#    migration) ──
quantities(::MagnetizationResponse) = (Magnetization,)
quantities(::SusceptibilityResponse) = (Susceptibility,)
quantities(::EntropyResponse) = (ThermalEntropy, FreeEnergy)
quantities(::GibbsHelmholtz) = (Energy, FreeEnergy)
quantities(::FreeEnergyFromZ) = (FreeEnergy, PartitionFunction)
quantities(::FreeEnergyLegendre) = (FreeEnergy, Energy, ThermalEntropy)
# Type-keyed relations (thermodynamic.jl, ensembles.jl, transport.jl): `quantities`
# is auto-derived from the typed slots, UNIONED with the `also_constrains` hints
# below. Those relations constrain a quantity through a SUPPLIED variance/derivative
# (not a typed subject), so the association is declared here to keep the physics
# graph complete (e.g. Var(M) → Magnetization in the susceptibility FDT):
also_constrains(::SusceptibilityFDT) = (Magnetization,)
also_constrains(::SpecificHeatFDT) = (Energy,)
also_constrains(::SpecificHeatFromEntropy) = (ThermalEntropy,)
also_constrains(::MottFormula) = (Conductivity,)   # transport: dlnσ/dε → Conductivity
also_constrains(::MicrocanonicalTemperature) = (ThermalEntropy, Energy)  # β = ∂S/∂E
# HeatCapacityDifference GAINS ThermalExpansionCoefficient + IsothermalCompressibility
# (α, κT are now typed subjects → auto). LinearResponseFDT + the 4 Maxwell relations
# stay symbol-keyed (generic / pure-derivative — no named quantity subject).

# ── Correlations / Green's functions ──
# NB: Dyson, SpectralFromGreens and the whole Keldysh block below are now
# TYPE-KEYED (`Name(x::Quantity, …)` in spectral.jl / keldysh.jl), so their
# `quantities` is auto-derived from the declaration — no hand-link here.  The
# remaining entries are legacy symbol-keyed relations awaiting migration.
quantities(::SpectralSumRule) = (SpectralFunction,)
quantities(::DetailedBalance) = (DynamicalStructureFactor,)
quantities(::DynamicalFDT) = (DynamicalStructureFactor, DynamicalSusceptibility)
quantities(::CorrelationLengthGap) = (CorrelationLength, MassGap)
quantities(::NMRExponent) = (NMRSpinRelaxationRate,)
quantities(::FiniteSizeGap) = (MassGap,)
function quantities(::StaticFromDynamicalStructureFactor)
    return (StaticStructureFactor, DynamicalStructureFactor)
end

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
quantities(::RenyiMonotonicity) = (RenyiEntropy,)
quantities(::RelativeEntropyNonNegativity) = (RelativeEntropy,)
quantities(::CFTEntanglementSlope) = (VonNeumannEntropy,)

# ── Quantum-mechanical foundations ──
# VirialTheorem is type-keyed (quantum.jl), `quantities` auto-derived. The Ehrenfest /
# Hellmann–Feynman / uncertainty / Lieb–Robinson relations stay symbol-keyed (generic
# operators / derivatives). EnergyVarianceEigenstate stays symbol-keyed — its subject
# EnergyVariance enters only through the moment combination ⟨H²⟩ − E², not a slot.
quantities(::EnergyVarianceEigenstate) = (EnergyVariance,)

# ── Topology ──
quantities(::TKNN) = (Conductivity,)
