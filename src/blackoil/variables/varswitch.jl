import Jutul: replace_value

function update_primary_variable!(state, pvar::BlackOilUnknown, state_symbol, model, Dx)
    v = state[state_symbol]
    rs_tab = model.system.saturation_table
    dr_max = 0.25
    # dr_max = pvar.dr_max
    ds_max = pvar.ds_max
    ϵ = 1e-6
    for i in eachindex(v)
        dx = Dx[i]
        old_x, old_state = v[i]
        if old_state == OilOnly
            p = state.Pressure[i]
            rs_sat = rs_tab(p)

            abs_rs_max = dr_max*rs_sat
            next_x = old_x + Jutul.choose_increment(value(old_x), dx, abs_rs_max, nothing, 0, nothing)
            if next_x > rs_sat
                # Switch to gas saturation as primary variable
                next_x = replace_value(next_x, ϵ)
                next_state = OilAndGas
            else
                next_state = old_state
            end
        else
            next_x = old_x + Jutul.choose_increment(value(old_x), dx, ds_max, nothing, nothing, 1)
            if next_x <= 0
                # Negative saturations - we switch to Rs as the primary variable

                p = state.Pressure[i]
                rs_sat = rs_tab(p)
                next_x = replace_value(next_x, rs_sat - ϵ)
                next_state = OilOnly
            else
                next_state = old_state
            end
        end
        v[i] = (next_x, next_state)
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
    return (x, state)
end

@jutul_secondary function update_as_secondary!(s, ph::BlackOilPhaseState, model::SimulationModel{D, S}, param, BlackOilUnknown)  where {D, S<:BlackOilVariableSwitchingSystem}
    for i in eachindex(s)
        s[i] = BlackOilUnknown[i][2]
    end
end

@jutul_secondary function update_as_secondary!(s, sat::Saturations, model::SimulationModel{D, S}, param, BlackOilUnknown, ImmiscibleSaturation) where {D, S<:BlackOilVariableSwitchingSystem}
    @assert size(s, 1) == 3
    T = eltype(s)
    for i in eachindex(BlackOilUnknown)
        sw = ImmiscibleSaturation[i]
        s[1, i] = sw
        x, phase_state = BlackOilUnknown[i]
        if phase_state == OilOnly
            sg = zero(T)
        else
            sg = x
        end
        s[2, i] = 1 - sw - sg
        s[3, i] = sg
    end
end


@jutul_secondary function update_as_secondary!(rs, ph::Rs, model::SimulationModel{D, S}, param, Pressure, BlackOilUnknown)  where {D, S<:BlackOilVariableSwitchingSystem}
    for i in eachindex(BlackOilUnknown)
        x, phase_state = BlackOilUnknown[i]
        if phase_state == OilOnly
            r = x
        else
            r = model.system.saturation_table(Pressure[i])
        end
        rs[i] = r
    end
end
