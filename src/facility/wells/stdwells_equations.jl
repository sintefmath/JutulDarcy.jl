function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::SimpleWellEquation, model, dt, ldisc = local_discretization(eq, self_cell)) where T_e
    # error()
    mass = state.MassFractions
    mass0 = state0.MassFractions
    p = state.Pressure[1]
    p0 = state0.Pressure[1]
    ϵ = 1e-20
    for i in eachindex(eq_buf)
        m = mass[i]
        m0 = mass0[i]
        if i == 1
            # Make sure that no sources/sinks means that the pressure remains
            # the same and mass fractions remain the same.
            Δ = abs(p - p0)
        else
            Δ = 0
        end
        eq_buf[i] = ϵ*(abs(m - m0) + Δ)
    end
end

function Jutul.convergence_criterion(model, storage, eq::SimpleWellEquation, eq_s, r; dt = 1)
    e = maximum(abs, r)
    R = (CNV = (errors = e, names = "R"), )
    return R
end
