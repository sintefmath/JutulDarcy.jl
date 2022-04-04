@enum PresentPhasesBlackOil OilOnly GasOnly OilAndGas


include("variables/variables.jl")
# include("utils.jl")
include("flux.jl")
# include("sources.jl")
include("wells.jl")

function select_primary_variables_system!(S, domain, system::BlackOilSystem, formulation)
    S[:Pressure] = Pressure(max_rel = 0.2, minimum = 1e5)
    S[:GasMassFraction] = GasMassFraction(dz_max = 0.1)
    if has_other_phase(system)
        S[:ImmiscibleSaturation] = ImmiscibleSaturation(ds_max = 0.2)
    end
end

function select_secondary_variables_system!(S, domain, system::BlackOilSystem, formulation)
    select_default_darcy!(S, domain, system, formulation)
    S[:Saturations] = Saturations()
    S[:PhaseState] = BlackOilPhaseState()
    S[:ShrinkageFactors] = ConstantVariables(1.0)
    S[:Rs] = Rs()
end

get_phases(sys::StandardBlackOilSystem{T, false}) where T = (LiquidPhase(), VaporPhase())
get_phases(sys::StandardBlackOilSystem) = (AqueousPhase(), LiquidPhase(), VaporPhase())
number_of_components(sys::StandardBlackOilSystem) = length(get_phases(sys))
