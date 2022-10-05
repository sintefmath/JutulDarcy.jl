export well_mismatch

function compute_well_qoi(well_model, state, well::Symbol, pos, rhoS, control)
    well_state = state[well]
    well_state = convert_to_immutable_storage(well_state)

    q_t = state[:Facility][:TotalSurfaceMassRate][pos]
    target = control.target

    rhoS, S = flash_wellstream_at_surface(well_model, well_state, rhoS)
    v = well_target(control, target, well_model, well_state, rhoS, S)
    if rate_weighted(target)
        v *= value(q_t)
    end
    return v
end

# function well_mismatch(model_f, model_c, states_f, states_c, tstep, forces, well::Symbol, qoi)
#     well_f = model_f[well]
#     well_c = model_c[well]
#     rhoS = reference_densities(well_f.system)
#     pos = get_well_position(model_f.models[:Facility].domain, well)

#     t_tot = sum(tstep)
#     mismatch = 0.0
#     for i in eachindex(tstep)
#         dt = tstep[i]
#         force = Jutul.forces_for_timestep(nothing, forces, tstep, i)
#         ctrl = force[:Facility].control[well]
#         ctrl = replace_target(ctrl, qoi)

#         state_f = states_f[i]
#         state_c = states_c[i]

#         qoi_f = compute_qoi(well_f, state_f, well, pos, rhoS, ctrl)
#         qoi_c = compute_qoi(well_c, state_c, well, pos, rhoS, ctrl)

#         mismatch += (dt/t_tot)*(qoi_f - qoi_c)^2
#     end
#     return sqrt(mismatch)
# end

function well_mismatch(qoi, wells, model_f, states_f, model_c, state_c, dt, step_no, forces; weights = ones(length(qoi)), scale = 1.0)
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

        state_f = states_f[step_no]

        for (i, q) in enumerate(qoi)
            ctrl = replace_target(ctrl, q)
            qoi_f = compute_well_qoi(well_f, state_f, well, pos, rhoS, ctrl)
            qoi_c = compute_well_qoi(well_c, state_c, well, pos, rhoS, ctrl)

            Δ = qoi_f - qoi_c
            obj += (weights[i]*Δ)^2
        end
    end
    return scale*dt*sqrt(obj)
end
