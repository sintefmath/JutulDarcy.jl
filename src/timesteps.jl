mutable struct ControlChangeTimestepSelector <: Jutul.AbstractTimestepSelector
    thresholds::Dict
    dt_after_change::Float64
    include_temperature::Bool
    prev_controls
    function ControlChangeTimestepSelector(thresholds, dt_after_change; include_temperature = false)
        new(thresholds, dt_after_change, include_temperature, nothing)
    end

end

function setup_control_change_timestep_selector(model, threshold = 0.25, dt_after_change = 5.0si_unit(:hour); wells = :all)

    all_wells = well_symbols(model)
    if wells == :all
        wells = all_wells
    else
        @assert all(well ∈ all_wells for well in wells) "Invalid well symbols"
    end
    if length(threshold) == 1
        threshold = fill(threshold, length(wells))
    end
    @assert length(threshold) == length(wells) 
        "Threshold must be provided either as a single value, or per well"

    thresholds = Dict{Symbol, Float64}()
    for (i, well) in enumerate(wells)
        @assert threshold[i] >= 0.0 "Threshold must be non-negative"
        thresholds[well] = threshold[i]
    end

    rmodel = reservoir_model(model)
    include_temperature = haskey(rmodel.primary_variables, :Temperature)

    ControlChangeTimestepSelector(
        thresholds, dt_after_change, include_temperature = include_temperature)

end

function Jutul.pick_next_timestep(sel::ControlChangeTimestepSelector, 
    sim, config, dt_prev, dT, forces, reports, current_reports, step_index, new_step)
    # Get the current controls
    curr_controls = forces[:Facility].control
    if isnothing(sel.prev_controls)
        sel.prev_controls = curr_controls
        return dT
    end
    change = false
    wells = keys(sel.thresholds)
    for well in wells
        # Check if the control type has changed
        ctrl0 = sel.prev_controls[well]
        ctrl = curr_controls[well]
        if typeof(ctrl) != typeof(ctrl0)
            change = true
            break
        end
        # If control is disabled, there are no changes to check
        if ctrl isa DisabledControl
            continue
        end
        # Check if the target type has changed
        target = ctrl.target
        target0 = ctrl0.target
        if typeof(ctrl.target) != typeof(ctrl0.target)
            change = true
            break
        end
        # Check if the target value has changed more than prescribed threshold
        threshold = sel.thresholds[well]
        val = target.value
        val0 = target0.value
        if abs(val - val0)/max(abs(val0), eps(1.0)) > threshold
            change = true
            break
        end
        # Check if the temperature has changed more than prescribed threshold
        if !(sel.include_temperature && ctrl isa InjectorControl)
            break
        end
        temp0 = sel.prev_controls[well].temperature
        temp = curr_controls[well].temperature
        if abs(temp - temp0)/max(temp0, eps(1.0)) > threshold
            change = true
            break
        end
    end
    # Update the previous controls
    sel.prev_controls = curr_controls
    # Update timetep if change has occured
    ΔT = change ? sel.dt_after_change : dT
    return ΔT
    
end