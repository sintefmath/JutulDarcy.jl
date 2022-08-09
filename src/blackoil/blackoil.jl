@enum PresentPhasesBlackOil OilOnly GasOnly OilAndGas


include("variables/variables.jl")
# include("utils.jl")
include("flux.jl")
# include("sources.jl")
include("wells.jl")

blackoil_formulation(::StandardBlackOilSystem{D, W, R, F}) where {D, W, R, F} = F

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
    select_default_darcy!(S, model.domain, model.system, model.formulation)
    S[:Saturations] = Saturations()
    S[:PhaseState] = BlackOilPhaseState()
    S[:ShrinkageFactors] = ConstantVariables(1.0)
    S[:Rs] = Rs()
end

get_phases(sys::StandardBlackOilSystem{T, false}) where T = (LiquidPhase(), VaporPhase())
get_phases(sys::StandardBlackOilSystem) = (AqueousPhase(), LiquidPhase(), VaporPhase())
number_of_components(sys::StandardBlackOilSystem) = length(get_phases(sys))
