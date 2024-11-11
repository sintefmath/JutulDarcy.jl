struct AMGXPreconditioner <: Jutul.JutulPreconditioner
    settings::Dict{String, Any}
    data::Dict{Symbol, Any}
    function AMGXPreconditioner(settings::Dict{String, Any})
        data = Dict{Symbol, Any}()
        new(settings, data)
    end
end

function AMGXPreconditioner(; kwarg...)
    settings = Dict{String, Any}()
    for (k, v) in kwarg
        settings["$k"] = v
    end
    return AMGXPreconditioner(settings)
end

const AMGXCPR = CPRPreconditioner{JutulDarcy.AMGXPreconditioner, <:Any}

function update_preconditioner!(amg::AMGXPreconditioner, A::Jutul.StaticSparsityMatrixCSR, b, context::ParallelCSRContext, executor)
    # Intentionally do nothing - a bit hackish
end

function JutulDarcy.gpu_update_preconditioner!(cpr::AMGXCPR, lsys, model, storage, recorder, executor, krylov, J_bsr, r_cu, op)
    @tic "CPU cpr work" Jutul.update_preconditioner!(cpr, lsys, model, storage, recorder, executor, update_system_precond = false)
    # Transfer pressure system to GPU
    @tic "update system precond" JutulDarcy.gpu_update_preconditioner!(cpr.system_precond, lsys, model, storage, recorder, executor, krylov, J_bsr, r_cu, op)
    @tic "update pressure system" JutulDarcy.update_amgx_pressure_system!(cpr.pressure_precond, cpr.storage.A_p, eltype(J_bsr), cpr, recorder)
    # How to get the linear operator in here?
    gpu_cpr_setup_buffers!(cpr, J_bsr, r_cu, op, recorder)
end

function update_amgx_pressure_system!

end

function gpu_cpr_setup_buffers!

end

function gpu_amgx_solve!

end

function Jutul.apply!(x, cpr::AMGXCPR, r)
    # Apply smoother
    @tic "system precond apply" Jutul.apply!(x, cpr.system_precond, r)
    # Correct the residual
    A_ps = cpr.pressure_precond.data[:operator]
    r_corrected = cpr.pressure_precond.data[:buffer_full]
    @tic "residual correction" JutulDarcy.correct_residual!(r_corrected, A_ps, x)
    # Construct pressure residual
    r_p = cpr.pressure_precond.data[:buffer_p]
    w_p = cpr.pressure_precond.data[:w_p]
    ncomp, ncell = size(w_p)
    @tic "residual reduction" gpu_reduce_residual!(r_p, w_p, r)
    # Apply pressure preconditioner
    amgx = cpr.pressure_precond
    @tic "AMGX apply" dp = gpu_amgx_solve!(amgx, r_p)
    # Update increments for pressure
    @tic "increment pressure" gpu_increment_pressure!(x, dp)
end
