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
                                    atol = nothing,
                                    v = 0,
                                    mode = :forward,
                                    solver = :bicgstab,
                                    max_iterations = nothing,
                                    update_interval = :iteration,
                                    update_interval_partial = :iteration,
                                    amg_type = default_amg_symbol(),
                                    smoother_type = :ilu0,
                                    max_coarse = 10,
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
    default_tol = 0.005
    max_it = 200
    if precond == :cpr
        if isnothing(cpr_type)
            if isa(model.system, ImmiscibleSystem)
                cpr_type = :analytical
            else
                cpr_type = :true_impes
            end
        end
        p_solve = default_psolve(max_coarse = max_coarse, type = amg_type)
        if smoother_type == :ilu0
            s = ILUZeroPreconditioner()
        elseif smoother_type == :jacobi
            s = JacobiPreconditioner()
        else
            error("Smoother :$smoother_type not supported for CPR.")
        end
        prec = CPRPreconditioner(p_solve, s, strategy = cpr_type, 
                                 update_interval = update_interval,
                                 partial_update = partial_update,
                                 update_interval_partial = update_interval_partial)
        default_tol = 1e-3
        max_it = 50
    elseif precond == :ilu0
        prec = ILUZeroPreconditioner()
    elseif precond == :jacobi
        prec = JacobiPreconditioner()
    elseif precond == :spai0
        prec = SPAI0Preconditioner()
    else
        error("Solver $precond not supported for $(model.context)")
    end
    if isnothing(rtol)
        rtol = default_tol
    end
    if isnothing(max_iterations)
        max_iterations = max_it
    end
    lsolve = GenericKrylov(solver, verbose = v, preconditioner = prec, 
            relative_tolerance = rtol, absolute_tolerance = atol,
            max_iterations = max_iterations; kwarg...)
    return lsolve
end

function default_amg_symbol()
    if Jutul.check_hypre_availability(throw = false)
        amg_type = :hypre
    else
        amg_type = :smoothed_aggregation
    end
    return amg_type
end


function Jutul.select_linear_solver(m::SimulationModel{<:Any, S, <:Any, <:Any}; kwarg...) where S<:MultiPhaseSystem
    return reservoir_linsolve(m; kwarg...)
end
