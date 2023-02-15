function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::SimpleWellEquation, model, dt, ldisc = local_discretization(eq, self_cell)) where T_e
    # error()
    mass = state.MassFractions
    mass0 = state0.MassFractions
    p = state.Pressure[1]
    p0 = state0.Pressure[1]
    ϵ = 1e-3
    for i in eachindex(eq_buf)
        m = mass[i]
        m0 = mass0[i]
        if i == 1
            Δ = (p - p0)*1e-18
        else
            Δ = 0
        end
        eq_buf[i] = ϵ*(m - m0)/dt + Δ
    end
end

function Jutul.convergence_criterion(model, storage, eq::SimpleWellEquation, eq_s, r; dt = 1)
    # For CNV = 0.001 we assume 1% Inf error in well composition balance
    ϵ = 0.01
    e = maximum(abs, r)
    scale = dt/ϵ
    R = (CNV = (errors = e*scale, names = "R"), )
    return R
end
