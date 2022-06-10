import Jutul: replace_value

function update_primary_variable!(state, pvar::BlackOilUnknown, state_symbol, model, Dx)
    v = state[state_symbol]
    rs_tab = model.system.saturation_table
    dr_max = pvar.dr_max
    ds_max = pvar.ds_max
    ϵ = 1e-4

    @inbounds for i in eachindex(v)
        dx = Dx[i]
        old_x, old_state, was_near_bubble = v[i]
        p = state.Pressure[i]
        rs_sat = rs_tab(value(p))
        if old_state == OilOnly
            abs_rs_max = dr_max*rs_sat
            next_x = old_x + Jutul.choose_increment(value(old_x), dx, abs_rs_max, nothing, 0, nothing)
            if next_x > rs_sat
                if was_near_bubble
                    # We are sufficiently close to the saturated point. Switch to gas saturation as primary variable.
                    # @info "$i Switching to saturated" value(next_x) value(old_x) value(rs_sat)
                    next_x = replace_value(next_x, ϵ)
                    next_state = OilAndGas
                    is_near_bubble = true
                else
                    # We are passing the saturated point, but we were sufficiently far from it that we limit the update
                    # to just before the saturated point.
                    next_x = replace_value(next_x, rs_sat*(1 - ϵ))
                    next_state = old_state
                    is_near_bubble = true
                end
            else
                next_state = old_state
                is_near_bubble = false
            end
        else
            next_x = old_x + Jutul.choose_increment(value(old_x), dx, ds_max, nothing, nothing, 1)
            if next_x <= 0
                if was_near_bubble
                    # Negative saturations - we switch to Rs as the primary variable
                    # @info "$i Switching to undersaturated" value(next_x) value(old_x)
                    next_x = replace_value(next_x, rs_sat*(1 - ϵ))
                    next_state = OilOnly
                    is_near_bubble = true
                else
                    next_x = replace_value(next_x, ϵ)
                    next_state = old_state
                    is_near_bubble = true
                end
            else
                next_state = old_state
                is_near_bubble = false
            end
        end
        # if next_state == OilOnly
        #     @assert 0 <= next_x <= rs_sat "$next_x $rs_sat"
        # else
        #     @assert 0 <= next_x <= 1 "$next_x"
        # end
        # if old_state != next_state
        #     @info "$i $old_state to $next_state" value(next_x) value(old_x) value(rs_sat) next_x
        # end
        v[i] = (next_x, next_state, is_near_bubble)
    end
end

function blackoil_unknown_init(F_rs, sg, rs, p)
    rs_sat = F_rs(p)
    if sg > 0
        @assert rs ≈ rs_sat
        x = sg
        state = OilAndGas
    else
        x = rs
        state = OilOnly
    end
    return (x, state, false)
end

@jutul_secondary function update_as_secondary!(s, ph::BlackOilPhaseState, model::SimulationModel{D, S}, param, BlackOilUnknown)  where {D, S<:BlackOilVariableSwitchingSystem}
    mb = minbatch(model.context)
    @inbounds @batch minbatch = mb for i in eachindex(s)
        s[i] = BlackOilUnknown[i][2]
    end
end

@jutul_secondary function update_as_secondary!(s, sat::Saturations, model::SimulationModel{D, S}, param, BlackOilUnknown, ImmiscibleSaturation) where {D, S<:BlackOilVariableSwitchingSystem}
    @assert size(s, 1) == 3
    T = eltype(s)
    mb = minbatch(model.context)
    @inbounds @batch minbatch = mb for i in eachindex(BlackOilUnknown)
        sw = ImmiscibleSaturation[i]
        s[1, i] = sw
        x, phase_state, = BlackOilUnknown[i]
        if phase_state == OilOnly
            sg = zero(T)
        else
            sg = x
        end
        s[2, i] = (1 - sw)*(1 - sg)
        s[3, i] = (1 - sw)*sg
    end
end


@jutul_secondary function update_as_secondary!(rs, ph::Rs, model::SimulationModel{D, S}, param, Pressure, BlackOilUnknown)  where {D, S<:BlackOilVariableSwitchingSystem}
    mb = minbatch(model.context)
    @inbounds @batch minbatch = mb for i in eachindex(BlackOilUnknown)
        x, phase_state, = BlackOilUnknown[i]
        if phase_state == OilOnly
            r = x
        else
            r = model.system.saturation_table(Pressure[i])
        end
        rs[i] = r
    end
end
