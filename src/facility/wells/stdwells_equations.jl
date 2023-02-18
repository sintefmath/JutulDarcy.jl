function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::SimpleWellEquation, model, dt, ldisc = local_discretization(eq, self_cell)) where T_e
    vol = model.domain.grid.volume
    mass = state.MassFractions
    mass0 = state0.MassFractions

    density(P) = 1.0 + (P - DEFAULT_MINIMUM_PRESSURE)*1e-8
    p = state.Pressure[1]
    p0 = state0.Pressure[1]

    rho = density(p)
    rho0 = density(p0)

    V = vol*density(p)
    V0 = vol*density(p0)
    for i in eachindex(eq_buf)
        m = mass[i]*V
        m0 = mass0[i]*V0
        eq_buf[i] = (m - m0)/dt
    end
end

function Jutul.convergence_criterion(model, storage, eq::SimpleWellEquation, eq_s, r; dt = 1)
    vol = model.domain.grid.volume
    scale = 0.1
    e = map(x -> scale*abs(x)*dt/vol, vec(r))
    R = (CNV = (errors = e, names = map(x -> "M$x", eachindex(e))), )
    return R
end
