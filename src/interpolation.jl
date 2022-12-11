import Jutul: first_lower, interval_weight, linear_interp

function bin_interval(t, x)
    ix = Jutul.first_lower(t, x)
    w = interval_weight(t, x, ix)
    return (w, ix)
end

function interp_pvt(pvto, p, v, type = :shrinkage; cap = false)
    P, R, SP, pos = pvt_table_vectors(pvto)
    pos = pvto.pos
    if cap
        v_max = linear_interp(SP, R, p)
        v = min(v_max, v)
    end
    w, ix = bin_interval(R, v)
    w = max(w, 0)
    # We now know what lines (for given v) bound the point
    # Get the positions of those lines in the linear array
    lower = pos[ix]:(pos[ix+1]-1)
    upper = pos[ix+1]:(pos[ix+2]-1)

    P_l = view(P, lower)
    P_u = view(P, upper)
    # Width of interval in saturation table
    @inbounds Δp = SP[ix+1] - SP[ix]

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

    return f_l*(1-w) + w*f_u
end
