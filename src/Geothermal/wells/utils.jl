using JutulDarcy

function JutulDarcy.InjectorControl(well_model, target::T, mix;
        enthalpy = missing,
        kwargs...
    ) where {T<:JutulDarcy.WellTarget}
    
    ismissing(enthalpy) || error("Enthaply cannot be defined for this constructor")

    function well_enthalpy(p, T)
        q = 0.0
        tp = Base.promote_type(eltype(p), eltype(T))

        state_tmp = setup_state(well_model, Pressure = p, Temperature = T, TotalMassFlux = q)
        for (k, v) in pairs(state_tmp)
            if k == :SurfaceWellConditions
                continue
            end
            if v isa Vector
                state_tmp[k] = Base.convert(Vector{tp}, v)
            elseif v isa AbstractMatrix
                state_tmp[k] = Base.convert(AbstractMatrix{tp}, v)
            else
                state_tmp[k] = Base.convert(tp, v)
            end
        end
        state_tmp = Jutul.evaluate_all_secondary_variables(well_model, state_tmp)
        h = state_tmp[:FluidEnthalpy][1]
        return h
    end

    return JutulDarcy.InjectorControl(target, mix; enthalpy = well_enthalpy, kwargs...)
    
end