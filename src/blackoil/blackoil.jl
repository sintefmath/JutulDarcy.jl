@enum PresentPhasesBlackOil OilOnly GasOnly OilAndGas


include("variables/variables.jl")
include("flux.jl")
include("wells.jl")
include("data.jl")

blackoil_formulation(::StandardBlackOilSystem{V, D, W, R, F}) where {V, D, W, R, F} = F

function select_primary_variables!(S, system::BlackOilSystem, model)
    S[:Pressure] = Pressure(max_rel = 0.2, minimum = 1e5)
    bf = blackoil_formulation(system)
    if bf == :varswitch
        S[:BlackOilUnknown] = BlackOilUnknown()
    elseif bf == :zg
        S[:GasMassFraction] = GasMassFraction(dz_max = 0.1)
    else
        error("Unsupported formulation $bf.")
    end
    if has_other_phase(system)
        S[:ImmiscibleSaturation] = ImmiscibleSaturation(ds_max = 0.2)
    end
end

function select_secondary_variables!(S, system::BlackOilSystem, model)
    select_default_darcy_secondary_variables!(S, model.domain, model.system, model.formulation)
    S[:Saturations] = Saturations()
    S[:PhaseState] = BlackOilPhaseState()
    spe1_data = blackoil_bench_pvt(:spe1)
    pvt = spe1_data[:pvt]
    S[:PhaseMassDensities] = DeckDensity(pvt)
    S[:ShrinkageFactors] = DeckShrinkageFactors(pvt)
    if !(model.domain.grid isa WellGrid)
        S[:SurfaceVolumeMobilities] = SurfaceVolumeMobilities()
    end
    S[:PhaseViscosities] = DeckViscosity(pvt)
    if has_disgas(system)
        S[:Rs] = Rs()
    end
    if has_vapoil(system)
        S[:Rv] = Rv()
    end
end

get_phases(sys::StandardBlackOilSystem) = sys.phases
number_of_components(sys::StandardBlackOilSystem) = length(get_phases(sys))
phase_indices(sys::StandardBlackOilSystem) = sys.phase_indices

has_vapoil(::Any) = false
has_disgas(::Any) = false

has_vapoil(::StandardBlackOilSystem) = true
has_disgas(::StandardBlackOilSystem) = true

has_vapoil(::DisgasBlackOilSystem) = false
has_disgas(::VapoilBlackOilSystem) = false
