"""
    reservoir_linsolve(model, precond = :cpr; <keyword arguments>)

Set up iterative linear solver for a reservoir model from [`setup_reservoir_model`](@ref).

# Arguments
- `model`: Reservoir model that will linearize the equations for the linear solver
- `precond=:cpr`: Preconditioner type to use: Either :cpr (Constrained-Pressure-Residual) or :ilu0 (block-incomplete-LU) (no effect if `solver = :direct`).
- `v=0`: verbosity (can lead to a large amount of output)
- `solver=:bicgstab`: the symbol of a Krylov.jl solver (typically :gmres or :bicgstab)
- `update_interval=:once`: how often the CPR AMG hierarchy is reconstructed (:once, :iteration, :ministep, :step)
- `update_interval_partial=:iteration`: how often the pressure system is updated in CPR
- `max_coarse`: max size of coarse level if using AMG
- `cpr_type=nothing`: type of CPR (`:true_impes`, `:quasi_impes` or `nothing` for automatic)
- `partial_update=true`: perform partial update of CPR preconditioner outside of AMG update (see above)
- `rtol=1e-3`: relative tolerance for the linear solver
- `max_iterations=100`: limit for linear solver iterations

Additional keywords are passed onto the linear solver constructor.
"""
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
    if solver == :lu
        return LUSolver()
    end
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

function Jutul.select_linear_solver(m::SimulationModel{<:Any, S, <:Any, <:Any}; kwarg...) where S<:MultiPhaseSystem
    return reservoir_linsolve(m; kwarg...)
end
