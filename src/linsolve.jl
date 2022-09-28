function reservoir_linsolve(model,  precond = :cpr;
                                    rtol = nothing,
                                    v = 0,
                                    solver = :bicgstab,
                                    update_interval = :once,
                                    update_interval_partial = :iteration,
                                    amg_type = :smoothed_aggregation,
                                    max_coarse = nothing,
                                    cpr_type = nothing,
                                    partial_update = update_interval == :once,
                                    kwarg...)
    model = reservoir_model(model)
    if !Jutul.is_cell_major(matrix_layout(model.context))
        return nothing
    end
    if precond == :cpr
        if isnothing(cpr_type)
            if isa(model.system, ImmiscibleSystem)
                cpr_type = :analytical
            else
                cpr_type = :true_impes
            end
        end
        if isnothing(max_coarse)
            max_coarse = Int64(ceil(0.05*number_of_cells(model.domain)))
            max_coarse = min(1000, max_coarse)
            max_coarse = 10
        end
        p_solve = default_psolve(max_coarse = max_coarse, type = amg_type)
        prec = CPRPreconditioner(p_solve, strategy = cpr_type, 
                                 update_interval = update_interval,
                                 partial_update = partial_update,
                                 update_interval_partial = update_interval_partial)
        if isnothing(rtol)
            rtol = 1e-3
        end
    elseif precond == :ilu0
        prec = ILUZeroPreconditioner()
        if isnothing(rtol)
            rtol = 0.005
        end
    else
        error("Solver $precond not supported for $(model.context)")
    end
    max_it = 200
    atol = 0.0

    lsolve = GenericKrylov(solver, verbose = v, preconditioner = prec, 
            relative_tolerance = rtol, absolute_tolerance = atol,
            max_iterations = max_it; kwarg...)
    return lsolve
end

