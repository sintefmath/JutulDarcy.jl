"""
    compute_well_qoi(model::MultiModel, state, forces, well::Symbol, target::Union{WellTarget, Type})

Compute the quantity of interest (QoI) for a specified well in a reservoir simulation.

# Arguments
- `model::MultiModel`: The simulation model (from `setup_reservoir_model`).
- `state`: The current state of the simulation.
- `forces`: The forces applied in the simulation.
- `well::Symbol`: The symbol representing the well for which the QoI is computed.
- `target::Union{WellTarget, Type}`: The target type or specific well target for the QoI computation.

# Returns
- The computed QoI for the specified well.

"""
function compute_well_qoi(model::MultiModel, state, forces, well::Symbol, target::Union{WellTarget, Type})
    well_model = model[well]
    rhoS = reference_densities(well_model.system)

    if haskey(model.models, :Facility)
        pos = get_well_position(model.models[:Facility].domain, well)
        ctrl = forces[:Facility].control[well]
    else
        pos = 1
        ctrl = forces[Symbol("$(well)_ctrl")].control[well]
    end

    translate_target_to_symbol
    if ctrl isa DisabledControl
        qoi = 0.0
    else
        if target isa Type
            if target<:SurfaceVolumeTarget
                if ctrl isa InjectorControl
                    tv = 1.0
                else
                    tv = -1.0
                end
            else
                tv = 100e5
            end
            target = target(tv)
        end
        ctrl = replace_target(ctrl, target)
        qoi = compute_well_qoi(well_model, state, well::Symbol, pos, rhoS, ctrl)
    end
    return qoi
end

function compute_well_qoi(well_model, state, well::Symbol, pos, rhoS, control)
    well_state = state[well]
    well_state = convert_to_immutable_storage(well_state)

    q_t = state[:Facility][:TotalSurfaceMassRate][pos]
    target = control.target

    rhoS, S = surface_density_and_volume_fractions(well_state)
    v = well_target(control, target, well_model, well_state, rhoS, S)
    if rate_weighted(target)
        v *= q_t
    end
    return v
end


"""
    well_mismatch(qoi, wells, model_f, states_f, model_c, state_c, dt, step_no, forces; <keyword arguments>)

Compute well mismatch for a set of qoi's (well targets) and a set of well symbols.
"""
function well_mismatch(qoi, wells, model_f, states_f, model_c, state_c, dt, step_no, forces; weights = ones(length(qoi)), scale = 1.0, signs = nothing)
    if !(qoi isa AbstractArray)
        qoi = [qoi]
    end
    if !(wells isa AbstractArray)
        wells = [wells]
    end
    obj = 0.0
    @assert length(weights) == length(qoi)
    for well in wells
        pos = get_well_position(model_c.models[:Facility].domain, well)

        well_f = model_f[well]
        well_c = model_c[well]
        rhoS = reference_densities(well_f.system)

        ctrl = forces[:Facility].control[well]
        if ctrl isa DisabledControl
            continue
        end

        state_f = states_f[step_no]

        for (i, q) in enumerate(qoi)
            ctrl = replace_target(ctrl, q)
            if !isnothing(signs)
                s = signs[i]
                if ctrl isa ProducerControl
                    sgn = -1
                else
                    sgn = 1
                end
                if s != sgn && s != 0
                    continue
                end
            end
            qoi_f = compute_well_qoi(well_f, state_f, well, pos, rhoS, ctrl)
            qoi_c = compute_well_qoi(well_c, state_c, well, pos, rhoS, ctrl)

            Δ = qoi_f - qoi_c
            obj += (weights[i]*Δ)^2
        end
    end
    return scale*dt*obj
end
