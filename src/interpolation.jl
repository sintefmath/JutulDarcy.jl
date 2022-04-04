
function bin_interval(t, x)
    ix = first_lower(t, x)
    w = interval_weight(t, x, ix)
    return (w, ix)
end

function interp_pvt(pvto, p, v, type = :shrinkage; cap = false)
    # v_tab = pvto.key
    v_tab = second_key(pvto)
    pos = pvto.pos
    if cap
        v_max = linear_interp(pvto.sat_pressure, v_tab, p)
        v = min(v_max, v)
    end
    w_l, ix = bin_interval(v_tab, v)
    # We now know what lines (for given v) bound the point

    lower = pos[ix]:pos[ix+1]-1
    upper = pos[ix+1]:pos[ix+2]-1

    P_l = view(pvto.pressure, lower)
    P_u = view(pvto.pressure, upper)
    u = first_lower(P_u, p)
    @inbounds Δp_u = P_u[u+1] - P_u[u]

    Δp = Δp_u
    w = w_l

    p_u = p + (1-w)*Δp
    p_l = p - w*Δp

    if type == :shrinkage
        F_tab = pvto.shrinkage
    else
        F_tab = pvto.viscosity
    end

    F_u = view(F_tab, upper)
    F_l = view(F_tab, lower)

    f_l = linear_interp(P_l, F_l, p_l)
    f_u = linear_interp(P_u, F_u, p_u)

    f = f_l*(1-w) + w*f_u
 
    return f
end
