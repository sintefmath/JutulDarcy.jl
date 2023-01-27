import Jutul: first_lower, interval_weight, linear_interp

function bin_interval(t, x)
    ix = Jutul.first_lower(t, x)
    w = interval_weight(t, x, ix)
    return (w, ix)
end

function interp_pvt(pvto, p, v, F_tab = pvt.shrinkage; cap = false)
    P, R, SP, pos = pvt_table_vectors(pvto)
    pos = pvto.pos
    if cap
        v_max = linear_interp(SP, R, p)
        v = min(v_max, v)
    end
    w, ix = bin_interval(R, v)
    # We now know what lines (for given v) bound the point
    # Get the positions of those lines in the linear array
    ix_start = @inbounds pos[ix]
    ix_middle = @inbounds pos[ix+1]
    ix_end = @inbounds pos[ix+2]
    lower = ix_start:(ix_middle-1)
    upper = ix_middle:(ix_end-1)

    P_l = @inbounds view(P, lower)
    P_u = @inbounds view(P, upper)
    # Width of interval in saturation table
    @inbounds Δp = SP[ix+1] - SP[ix]

    p_u = p + (1-w)*Δp
    p_l = p - w*Δp

    F_u = @inbounds view(F_tab, upper)
    F_l = @inbounds view(F_tab, lower)

    f_l = linear_interp(P_l, F_l, p_l)
    f_u = linear_interp(P_u, F_u, p_u)

    return f_l*(1.0-w) + w*f_u
end
