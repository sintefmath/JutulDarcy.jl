const MSWellDomain = DiscretizedDomain{<:MultiSegmentWell}
const MSWellFlowModel = SimulationModel{<:MSWellDomain, <:MultiPhaseSystem}

function Jutul.initialize_extra_state_fields!(state, domain::MSWellDomain, model; T = Float64)
    state[:SegmentModels] = copy(physical_representation(domain).segment_models)
end

function Jutul.update_before_step!(storage, domain::MSWellDomain, model, dt, forces; kwarg...)
    defaults = physical_representation(domain).segment_models
    models = isnothing(forces) ? defaults : something(get(forces, :segment_models, nothing), defaults)
    copyto!(storage.state.SegmentModels, models)
    return nothing
end

# Selection of primary variables
function select_primary_variables!(S, domain::MSWellDomain, model::MSWellFlowModel)
    if count_active_entities(domain, Faces()) > 0
        S[:TotalMassFlux] = TotalMassFlux()
    end
end

function select_equations!(eqs, domain::MSWellDomain, model::MSWellFlowModel)
    if count_active_entities(domain, Faces()) > 0
        w = physical_representation(domain)
        eqs[:potential_balance] = PotentialDropBalanceWell(
            domain.discretizations.mass_flow; include_inertia = w.include_inertia)
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
    return vars
end

function get_neighborship(W::MultiSegmentWell)
    return W.neighborship
end
