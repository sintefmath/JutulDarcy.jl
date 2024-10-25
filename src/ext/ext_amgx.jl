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

function update_preconditioner!(amg::AMGXPreconditioner, A::Jutul.StaticSparsityMatrixCSR, b, context::ParallelCSRContext, executor)
    # Intentionally do nothing - a bit hackish
end

function JutulDarcy.gpu_update_preconditioner!(cpr::CPRPreconditioner{JutulDarcy.AMGXPreconditioner, <:Any}, lsys, model, storage, recorder, executor, krylov, J_bsr, r_cu)
    Jutul.update_preconditioner!(cpr, lsys, model, storage, recorder, executor)
    
    JutulDarcy.update_amgx_pressure_system!(cpr.pressure_precond, cpr.storage.A_p)
    # Next transfer the coarse system to GPU
    # Update if needed
    # Then do the whole apply here
    error()
    # Jutul.update_preconditioner!(krylov.preconditioner, J_bsr, r_cu, model.context, executor)
end

function update_amgx_pressure_system!

end