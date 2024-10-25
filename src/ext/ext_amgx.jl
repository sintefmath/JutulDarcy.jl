struct AMGXPreconditioner <: Jutul.JutulPreconditioner
    settings::Dict{Symbol, Any}
    data::Dict{Symbol, Any}
    function AMGXPreconditioner(settings::Dict{Symbol, Any})
        data = Dict{Symbol, Any}()
        new(settings, data)
    end
end

function AMGXPreconditioner(; kwarg...)
    settings = Dict{Symbol, Any}()
    for (k, v) in kwarg
        settings[k] = v
    end
    return AMGXPreconditioner(settings)
end

const AMGXCPR = CPRPreconditioner{JutulDarcy.AMGXPreconditioner, <:Any}

function update_preconditioner!(amg::AMGXPreconditioner, A::Jutul.StaticSparsityMatrixCSR, b, context::ParallelCSRContext, executor)
    # Intentionally do nothing - a bit hackish
end

function JutulDarcy.gpu_update_preconditioner!(cpr::AMGXCPR, lsys, model, storage, recorder, executor, krylov, J_bsr, r_cu)
    Jutul.update_preconditioner!(cpr, lsys, model, storage, recorder, executor, update_system_precond = false)
    # Transfer pressure system to GPU
    JutulDarcy.update_amgx_pressure_system!(cpr.pressure_precond, cpr.storage.A_p)
end

function update_amgx_pressure_system!

end

function Jutul.apply!(x, cpr::AMGXCPR, r)
    # Apply smoother
    cpr_s = cpr.storage
    bz = cpr_s.block_size

    Jutul.apply!(x, cpr.system_precond, r)
    # Correct the residual
    # Construct pressure residual
    # Apply pressure preconditioner
    # Update increments for pressure
    error()
end
