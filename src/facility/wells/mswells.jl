const MSWellDomain = DiscretizedDomain{<:MultiSegmentWell}
const MSWellFlowModel = SimulationModel{<:MSWellDomain, <:MultiPhaseSystem}
# Selection of primary variables
function select_primary_variables!(S, ::MSWellDomain, model::MSWellFlowModel)
    S[:TotalMassFlux] = TotalMassFlux()
end

function select_equations!(eqs, domain::MSWellDomain, model::MSWellFlowModel)
    eqs[:potential_balance] = PotentialDropBalanceWell(domain.discretizations.mass_flow)
end

function select_parameters!(prm, domain::MSWellDomain, model::MSWellFlowModel)
    prm[:WellIndices] = WellIndices()
    prm[:PerforationGravityDifference] = PerforationGravityDifference()
    # TODO: Only add these params if model has thermal
    prm[:WellIndicesThermal] = WellIndicesThermal()
    prm[:MaterialThermalConductivities] = MaterialThermalConductivities()
    prm[:MaterialHeatCapacities] = MaterialHeatCapacities()
    prm[:MaterialDensities] = MaterialDensities()
end

function Jutul.select_minimum_output_variables!(vars, domain::DiscretizedDomain, model::MSWellFlowModel)
    push!(vars, :PhaseMassDensities)
    push!(vars, :Saturations)
    push!(vars, :SurfaceWellConditions)
    return vars
end

function get_neighborship(W::MultiSegmentWell)
    return W.neighborship
end
