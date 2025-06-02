function default_value(model, ::SurfaceWellConditions)
    rho = reference_densities(flow_system(model.system))
    return TopConditions(length(rho), density = rho, volume_fractions = missing)
end

function initialize_variable_value(model, pvar::SurfaceWellConditions, val::AbstractDict; need_value = false, T = Jutul.float_type(model.context))
    @assert need_value == false
    initialize_variable_value(model, pvar, [default_value(model, pvar)], T = T)
end

function initialize_variable_value(model, pvar::SurfaceWellConditions, val::Vector; need_value = false, T = Jutul.float_type(model.context))
    @assert need_value == false
    if length(val) > 1 
        @warn "Expected a single value, got $(length(val))"
    end
    tc = val[1]
    if eltype(tc.density) != T
        density = T.(tc.density)
        volume_fractions = T.(tc.volume_fractions)
        tc = TopConditions(density, volume_fractions)
    end
    return [tc]
end

function update_secondary_variable!(x::Vector{TopConditions{N, R}}, var::SurfaceWellConditions, model, state, ix) where {N, R}
    nstages = length(var.separator_conditions)
    rhoS = vol = missing
    if nstages == 0
        rhoS = reference_densities(model.system)
        cond = physical_representation(model).surface
        rhoS, vol = flash_wellstream_at_surface(var, model, model.system, state, rhoS, cond)
    else
        rhoS, vol = separator_surface_flash!(var, model, model.system, state)
    end
    x[1] = TopConditions(N, R, density = rhoS, volume_fractions = vol)
end

function Jutul.default_values(model, var::SurfaceWellConditions)
    return [default_value(model, var)]
end

function Jutul.get_dependencies(x::SurfaceWellConditions, model)
    return [:TotalMasses]
end

function initialize_variable_ad!(state, model, pvar::SurfaceWellConditions, symb, npartials, diag_pos; context = DefaultContext(), kwarg...)
    v_ad = get_ad_entity_scalar(1.0, npartials, diag_pos; kwarg...)
    ∂T = typeof(v_ad)
    nph = number_of_phases(model.system)
    state[symb] = [TopConditions(nph, ∂T)]
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
    vals
end

function add_separator_stage!(var::SurfaceWellConditions, cond = default_surface_cond(), dest = (0, 0); clear = false)
    sc = var.separator_conditions
    t = var.separator_targets

    @assert cond.p > 0.0
    @assert cond.T > 0.0
    for i in dest
        @assert i >= 0
    end
    if clear
        empty!(sc)
        empty!(t)
    end
    @assert length(t) == length(sc)
    push!(sc, cond)
    push!(t, dest)
    return var
end

function add_separator_stage!(model::SimulationModel, arg...; kwarg...)
    add_separator_stage!(model[:SurfaceWellConditions], arg...; kwarg...)
    return model
end
