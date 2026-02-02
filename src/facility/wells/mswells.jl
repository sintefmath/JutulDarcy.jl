const MSWellDomain = DiscretizedDomain{<:MultiSegmentWell}
const MSWellFlowModel = SimulationModel{<:MSWellDomain, <:MultiPhaseSystem}
# Selection of primary variables
function select_primary_variables!(S, domain::MSWellDomain, model::MSWellFlowModel)
    if count_active_entities(domain, Faces()) > 0
        S[:TotalMassFlux] = TotalMassFlux()
    end
end

function select_equations!(eqs, domain::MSWellDomain, model::MSWellFlowModel)
    if count_active_entities(domain, Faces()) > 0
        eqs[:potential_balance] = PotentialDropBalanceWell(domain.discretizations.mass_flow)
    end
end

function select_parameters!(prm, domain::MSWellDomain, model::MSWellFlowModel)
    prm[:WellIndices] = WellIndices()
    prm[:PerforationGravityDifference] = PerforationGravityDifference()
    if count_active_entities(domain, Faces()) > 0
        prm[:SegmentConnectionGravityDifference] = SegmentConnectionGravityDifference()
        prm[:SegmentRadius] = SegmentRadius()
        prm[:SegmentRadiusInner] = SegmentRadiusInner()
        prm[:SegmentCasingThickness] = SegmentCasingThickness()
        prm[:SegmentRoughness] = SegmentRoughness()
        prm[:SegmentLength] = SegmentLength()
    end
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
