function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::SimpleWellEquation, model, dt, ldisc = local_discretization(eq, self_cell)) where T_e
    # error()
    # mass = state.MassFractions
    ncomp = number_of_components(model.system)
    for i in 1:ncomp
        eq_buf[i] = 0.0
    end
end

function Jutul.convergence_criterion(model, storage, eq::SimpleWellEquation, eq_s, r; dt = 1)
    R = (CNV = (errors = maximum(abs, r)/dt, names = "R"), )
    return R
end
