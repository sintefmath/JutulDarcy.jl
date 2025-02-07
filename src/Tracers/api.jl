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

"""
    tracer_total_mass(tracer::AbstractTracer, model, state, concentration, resident_mass, cell, index)

# Arguments:
    - `tracer::AbstractTracer`: Tracer in question
    - `model`: Model object.
    - `state`: Current state
    - `concentration`: Tracer concentration.
    - `resident_mass_density`: Resident mass density (i.e. total mass density of all associated phases where the tracer is present)
    - `volume`: Fluid volume in cell.
    - `cell`: Cell index.
    - `index`: Index of tracer in the list of tracers
"""
function tracer_total_mass(tracer::AbstractTracer, model, state, concentration, resident_mass_density, volume, cell, index)
    return concentration*resident_mass_density*volume
end

function tracer_total_mass_outer(tracer::AbstractTracer, model, state, concentration, resident_mass_density, vol, cell, index)
    if resident_mass_density <= TRACER_TOL
        v = vol*Jutul.replace_value(concentration, 0.0)
    else
        v = tracer_total_mass(tracer, model, state, concentration, resident_mass_density, vol, cell, index)
    end
    return v
end

"""
    tracer_total_mass_flux(tracer, model, state, phase_mass_fluxes, C, index, upw, T = eltype(phase_mass_fluxes))

# Arguments:
    - `tracer`: Tracer in question
    - `model`: Model object.
    - `state`: Current state
    - `phase_mass_fluxes`: Mass fluxes for each phase
    - `C`: Tracer concentrations (as a matrix with the tracers in each cell as the columns)
    - `index`: Index of tracer in the list of tracers
    - `upw`: Upwind scheme
    - `T`: Type of the tracer concentrations
"""
function tracer_total_mass_flux(tracer, model, state, phase_mass_fluxes, index, upw, C = state.TracerConcentrations, T = eltype(C))
    v = zero(T)
    for phase in tracer_phase_indices(tracer)
        q_ph = phase_mass_fluxes[phase]
        C_iface = phase_upwind(upw, C, index, q_ph)
        v += C_iface*q_ph
    end
    return v
end

function add_tracers_to_model!(model::MultiModel, tracers; names = missing, kwarg...)
    if !isnothing(tracers)
        tracers = TracerFluxType(tracers, names = names)
        rmodel = reservoir_model(model)
        add_tracers_to_model!(rmodel, tracers; kwarg...)
        for (k, m) in pairs(model.models)
            if JutulDarcy.model_or_domain_is_well(m)
                # TODO: Add support for tracer for multisegment wells by
                # specializing the equation update.
                physical_representation(m.domain)::JutulDarcy.SimpleWell
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
