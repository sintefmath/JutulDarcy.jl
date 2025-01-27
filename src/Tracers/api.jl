function number_of_tracers(model::MultiModel)
    return number_of_tracers(reservoir_model(model))
end

function number_of_tracers(model::SimulationModel)
    t = get(model.equations, :tracers, nothing)
    if isnothing(t)
        n = 0
    else
        n = number_of_tracers(t.flux_type)
    end
    return n
end

number_of_tracers(t::TracerFluxType) = length(t.tracers)


function add_tracers_to_model!(model::MultiModel, tracers; names = missing, kwarg...)
    if !isnothing(tracers)
        tracers = TracerFluxType(tracers, names = names)
        rmodel = reservoir_model(model)
        add_tracers_to_model!(rmodel, tracers; kwarg...)
        for (k, m) in pairs(model.models)
            if JutulDarcy.model_or_domain_is_well(m)
                add_tracers_to_model!(m, tracers; kwarg...)
            end
        end
        # Add cross-terms here.
        new_ctp = []
        function add_next!(ctp, ctt)
            if ctp.target_equation == :mass_conservation
                te = :tracers
            else
                te = ctp.target_equation
            end
            if ctp.source_equation == :mass_conservation
                se = :tracers
            else
                se = ctp.source_equation
            end
            v = Jutul.CrossTermPair(ctp.target, ctp.source, te, se, ctt)
            push!(new_ctp, v)
        end
        for ctp in model.cross_terms
            ct = ctp.cross_term
            if ctp.cross_term isa JutulDarcy.ReservoirFromWellFlowCT
                ctt = ReservoirFromWellTracerCT(ct.WI, ct.reservoir_cells, ct.well_cells)
                add_next!(ctp, ctt)
            end
            if ct isa JutulDarcy.WellFromFacilityFlowCT
                add_next!(ctp, WellFromFacilityTracerCT(ct.well))
            end
        end
        for ctp in new_ctp
            add_cross_term!(model.cross_terms, ctp)
        end
    end
    return model
end

function add_tracers_to_model!(model::SimulationModel, tracers; names = missing, discretization = missing)
    if !isnothing(tracers)
        if ismissing(discretization)
            discretization = model.domain.discretizations.mass_flow
        end
        tracers = TracerFluxType(tracers, names = names)
        N = number_of_tracers(tracers)
        !haskey(model.equations, :tracers) || error("Tracers already set in model")
        model.equations[:tracers] = ConservationLaw(discretization, :TracerMasses, N, flux = tracers)
        model.secondary_variables[:TracerMasses] = TracerMasses(tracers)
        model.primary_variables[:TracerConcentrations] = TracerConcentrations(tracers)

        push!(model.output_variables, :TracerMasses)
        push!(model.output_variables, :TracerConcentrations)
        unique!(model.output_variables)
    end
    return model
end
