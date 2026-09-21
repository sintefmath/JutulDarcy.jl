function default_value(model::FacilityModel, ::SurfaceWellConditions)
    rho = reference_densities(model.system.multiphase)
    return TopConditions(length(rho), density = rho, volume_fractions = missing)
end

function initialize_variable_value(model::FacilityModel, pvar::SurfaceWellConditions, val::AbstractDict; need_value = false, T = Jutul.float_type(model.context))
    @assert need_value == false
    return initialize_variable_value(model, pvar, Jutul.default_values(model, pvar), T = T)
end

function initialize_variable_value(model::FacilityModel, pvar::SurfaceWellConditions, val::Vector; need_value = false, T = Jutul.float_type(model.context))
    @assert need_value == false
    n = number_of_entities(model, pvar)
    if length(val) == 1 && n > 1
        val = fill(only(val), n)
    else
        @assert length(val) == n "Expected $n surface-condition values, got $(length(val))"
    end
    return map(val) do tc
        if Jutul.numerical_type(typeof(tc)) == T
            tc
        else
            TopConditions(T.(tc.density), T.(tc.volume_fractions))
        end
    end
end

function update_secondary_variable!(x::Vector{TopConditions{N, R}}, var::SurfaceWellConditions, model::FacilityModel, state, ix) where {N, R}
    system = model.system.multiphase
    for well in ix
        rhoS, vol = surface_density_and_volume_fractions(var, model, system, state, well)
        x[well] = TopConditions(N, R, density = rhoS, volume_fractions = vol)
    end
    return x
end

function surface_density_and_volume_fractions(var, model, system::MultiPhaseSystem, state, well)
    rhoS = reference_densities(system)
    rates = @view state.SurfaceComponentRates[:, well]
    vol = rates./rhoS
    total_volume = sum(vol)
    if abs(value(total_volume)) < MIN_ACTIVE_WELL_RATE
        vol = fill(one(eltype(vol))/length(vol), length(vol))
    else
        vol = vol./total_volume
    end
    return (rhoS, vol)
end

function Jutul.default_values(model::FacilityModel, var::SurfaceWellConditions)
    return [default_value(model, var) for _ in 1:number_of_entities(model, var)]
end

function Jutul.get_dependencies(x::SurfaceWellConditions, model::FacilityModel)
    return (:SurfaceComponentRates,)
end

function initialize_variable_ad!(state, model::FacilityModel, pvar::SurfaceWellConditions, symb, npartials, diag_pos; context = DefaultContext(), kwarg...)
    v_ad = get_ad_entity_scalar(1.0, npartials, diag_pos; kwarg...)
    ∂T = typeof(v_ad)
    nph = number_of_phases(model.system.multiphase)
    nw = number_of_entities(model, pvar)
    state[symb] = [TopConditions(nph, ∂T) for _ in 1:nw]
    return state
end

function Jutul.numerical_type(::Type{TopConditions{N, T}}) where {N, T}
    return T
end

function Base.convert(::Type{TopConditions{N, Float64}}, v::TopConditions{N, <:ForwardDiff.Dual}) where N
    rho = value.(v.density)
    s = value.(v.volume_fractions)
    return TopConditions(N, Float64, density = rho, volume_fractions = s)
end

function Jutul.value(tc::TopConditions{N, <:ForwardDiff.Dual}) where N
    d = value.(tc.density)
    v = value.(tc.volume_fractions)
    return TopConditions(d, v)
end

@inline function Jutul.update_values!(vals::Vector{TopConditions{N, T}}, next::Vector{TopConditions{N, Float64}}) where {N, T<:ForwardDiff.Dual}
    for i in eachindex(vals)
        v0 = vals[i]
        v = next[i]
        rho = v0.density - value(v.density) + v.density
        vol = v0.volume_fractions - value(v.volume_fractions) + v.volume_fractions
        vals[i] = TopConditions(N, T, density = rho, volume_fractions = vol)
    end
    return vals
end

function add_separator_stage!(var::SurfaceWellConditions, cond = default_surface_cond(), dest = (0, 0); well = 1, clear = false)
    sc = var.separator_conditions[well]
    t = var.separator_targets[well]

    @assert cond.p > 0.0
    @assert cond.T > 0.0
    for i in dest
        @assert i >= 0
    end
    if clear
        empty!(sc)
        empty!(t)
        empty!(var.storage[well])
    end
    @assert length(t) == length(sc)
    push!(sc, cond)
    push!(t, dest)
    return var
end

function add_separator_stage!(model::FacilityModel, cond = default_surface_cond(), dest = (0, 0); well = only(model.domain.well_symbols), clear = false)
    pos = well isa Symbol ? get_well_position(model.domain, well) : well
    isnothing(pos) && throw(ArgumentError("Well $well is not controlled by this facility"))
    add_separator_stage!(model[:SurfaceWellConditions], cond, dest; well = pos, clear = clear)
    return model
end
