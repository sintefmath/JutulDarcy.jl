import Jutul: replace_value

function update_primary_variable!(state, pvar::BlackOilUnknown, state_symbol, model, Dx)
    v = state[state_symbol]
    sys = model.system
    rs_tab = sys.rs_max
    rv_tab = sys.rv_max
    dr_max = pvar.dr_max
    ds_max = pvar.ds_max
    keep_bub = sys.keep_bubble_flag
    sat_chop = sys.saturated_chop
    pressure = state.Pressure
    active = active_entities(model.domain, associated_entity(pvar), for_variables = true)
    update_bo_internal!(v, Dx, dr_max, ds_max, rs_tab, rv_tab, keep_bub, sat_chop, pressure, active)
end


function update_bo_internal!(v, Dx, dr_max, ds_max, rs_tab, rv_tab, keep_bub, sat_chop, pressure, active)
    @inbounds for (i, dx) in zip(active, Dx)
        varswitch_update_inner!(v, i, dx, dr_max, ds_max, rs_tab, rv_tab, keep_bub, sat_chop, pressure)
    end
end

Base.@propagate_inbounds function varswitch_update_inner!(v, i, dx, dr_max, ds_max, rs_tab, rv_tab, keep_bubble, sat_chop, pressure)
    ϵ = 1e-8
    X = v[i]
    old_x = X.val
    old_state = X.phases_present
    was_near_bubble = sat_chop || X.sat_close
    if old_state == OilAndGas
        next_x = old_x + Jutul.choose_increment(value(old_x), dx, ds_max)
        if next_x <= 0
            maybe_state = OilOnly
            tab = rs_tab
        elseif next_x >= 1
            maybe_state = GasOnly
            tab = rv_tab
        else
            maybe_state = OilAndGas
            tab = nothing
        end
        next_x, next_state, is_near_bubble = handle_phase_disappearance(pressure, i, tab, next_x, old_state, maybe_state, was_near_bubble, keep_bubble, ϵ)
    else
        if old_state == GasOnly
            tab = rv_tab
        else
            tab = rs_tab
        end
        next_x, next_state, is_near_bubble = handle_phase_appearance(pressure, i, tab, dr_max, old_state, old_x, dx, was_near_bubble, ϵ, keep_bubble)
    end
    v[i] = BlackOilX(next_x, next_state, is_near_bubble)
end

function handle_phase_disappearance(pressure, i, ::Nothing, next_x, old_state, possible_new_state, was_near_bubble, keep_bubble, ϵ)
    return (next_x, old_state, false)
end

function handle_phase_disappearance(pressure, i, r_tab, next_x, old_state, possible_new_state, was_near_bubble, keep_bubble, ϵ)
    if was_near_bubble
        # Gas/Oil Phase disappears and Rs/Rv becomes primary variable
        p = pressure[i]
        r_sat = r_tab(value(p))
        next_x = replace_value(next_x, r_sat*(1 - ϵ))
        next_state = possible_new_state
        is_near_bubble = keep_bubble
    else
        if possible_new_state == GasOnly
            next_sg = 1 - ϵ
        else
            next_sg = ϵ
        end
        next_x = replace_value(next_x, next_sg)
        next_state = old_state
        is_near_bubble = true
    end
    return (next_x, next_state, is_near_bubble)
end

function handle_phase_appearance(pressure, i, r_tab, dr_max, old_state, old_x, dx, was_near_bubble, ϵ, keep_bubble)
    p = pressure[i]
    r_sat = r_tab(value(p))

    abs_r_max = dr_max*r_sat
    next_x = old_x + Jutul.choose_increment(value(old_x), dx, abs_r_max, nothing, 0, nothing)
    if next_x > r_sat
        if was_near_bubble
            # We are sufficiently close to the saturated point. Switch to gas saturation as primary variable.
            if old_state == OilOnly
                sg = ϵ
            else
                sg = 1 - ϵ
            end
            next_x = replace_value(next_x, sg)
            next_state = OilAndGas
            is_near_bubble = keep_bubble
        else
            # We are passing the saturated point, but we were sufficiently far from it that we limit the update
            # to just before the saturated point.
            next_x = replace_value(next_x, r_sat*(1 - ϵ))
            next_state = old_state
            is_near_bubble = true
        end
    else
        next_state = old_state
        is_near_bubble = false
    end
    return (next_x, next_state, is_near_bubble)
end

function s_removed(s, d)
    if d < 1e-12
        s_bar = 1.0
    else
        s_bar = s/(1-d)
    end
    return s_bar
end

function blackoil_unknown_init(F_rs, F_rv::Nothing, sw, so, sg, rs, rv, p)
    rs_sat = F_rs(p)
    if sg > 0
        @assert rs ≈ rs_sat "rs = $rs is different from rs_sat = $rs_sat: sg = $sg so = $so"
        x = s_removed(sg, sw)
        state = OilAndGas
    else
        x = rs
        state = OilOnly
    end
    return BlackOilX(x, state, false)
end

function blackoil_unknown_init(F_rs::Nothing, F_rv, sw, so, sg, rs, rv, p)
    rv_sat = F_rv(p)
    if so > 0
        @assert rv ≈ rv_sat "rv = $rv is different from rv_sat = $rv_sat: sg = $sg so = $so"
        x = s_removed(sg, sw)
        state = OilAndGas
    else
        x = rv
        state = GasOnly
    end
    return BlackOilX(x, state, false)
end

function blackoil_unknown_init(F_rs, F_rv, sw, so, sg, rs, rv, p)
    rs_sat = F_rs(p)
    rv_sat = F_rv(p)
    if sg > 0 && so > 0
        @assert rv ≈ rv_sat "rv = $rv is different from rv_sat = $rv_sat: sg = $sg so = $so"
        @assert rs ≈ rs_sat "rs = $rs is different from rs_sat = $rs_sat: sg = $sg so = $so"
        x = s_removed(sg, sw)
        state = OilAndGas
    elseif so > 0
        x = rs
        state = OilOnly
    else
        x = rv
        state = GasOnly
    end
    return BlackOilX(x, state, false)
end

Jutul.value(t::BlackOilX{T}) where T<:ForwardDiff.Dual = BlackOilX(value(t.val), t.phases_present, t.sat_close)

function Jutul.update_values!(old::AbstractVector{BlackOilX}, new::AbstractVector{BlackOilX})
    for (i, v) in enumerate(new)
        o = old[i]
        old[i] = BlackOilX(o.val - value(o.val) + value(v.val), v.phases_present, v.sat_close)
    end
end

# Overloads for our specific data type
import Base:+
function (+)(l::BlackOilX, r::BlackOilX)
    return BlackOilX(l.val + r.val, r.phases_present, r.sat_close)
end

import Base:-
function (-)(l::BlackOilX, r::BlackOilX)
    return BlackOilX(l.val - r.val, r.phases_present, r.sat_close)
end

@jutul_secondary function update_phase_state!(s, ph::BlackOilPhaseState, model::SimulationModel{D, S}, BlackOilUnknown, ix)  where {D, S<:BlackOilVariableSwitchingSystem}
    @inbounds for i in ix
        s[i] = BlackOilUnknown[i].phases_present
    end
end

@jutul_secondary function updaet_saturations!(s, sat::Saturations, model::SimulationModel{D, S}, BlackOilUnknown, ImmiscibleSaturation, ix) where {D, S<:BlackOilVariableSwitchingSystem}
    @assert size(s, 1) == 3
    T = eltype(s)
    @inbounds for i in ix
        sw = ImmiscibleSaturation[i]
        s[1, i] = sw
        X = BlackOilUnknown[i]
        if X.phases_present == OilOnly
            sg = zero(T)
        elseif X.phases_present == GasOnly
            sg = one(T)
        else
            sg = X.val
        end
        s[2, i] = (1 - sw)*(1 - sg)
        s[3, i] = (1 - sw)*sg
    end
end


@jutul_secondary function update_rs!(rs, ph::Rs, model::SimulationModel{D, S}, Pressure, BlackOilUnknown, ix)  where {D, S<:BlackOilVariableSwitchingSystem}
    @inbounds for i in ix
        X = BlackOilUnknown[i]
        if X.phases_present == OilOnly
            r = X.val
        else
            p = @inbounds Pressure[i]
            r = model.system.rs_max(p)
        end
        rs[i] = r
    end
end

@jutul_secondary function update_rv!(rv, ph::Rv, model::SimulationModel{D, S}, Pressure, BlackOilUnknown, ix)  where {D, S<:BlackOilVariableSwitchingSystem}
    @inbounds for i in ix
        X = BlackOilUnknown[i]
        if X.phases_present == GasOnly
            r = X.val
        else
            p = @inbounds Pressure[i]
            r = model.system.rv_max(p)
        end
        rv[i] = r
    end
end
