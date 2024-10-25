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
