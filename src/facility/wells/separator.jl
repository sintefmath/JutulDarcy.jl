function default_value(model, ::SurfaceWellConditions)
    rho = reference_densities(model.system)
    return TopConditions(length(rho), density = rho, volume_fractions = missing)
end

function initialize_variable_value(model, pvar::SurfaceWellConditions, val::AbstractDict; need_value = false)
    @assert need_value == false
    initialize_variable_value(model, pvar, [default_value(model, pvar)])
end

function initialize_variable_value(model, pvar::SurfaceWellConditions, val::Vector; need_value = false)
    @assert need_value == false
    @assert length(val) == 1
    return val
end

function update_secondary_variable!(x::Vector{TopConditions{N, R}}, var::SurfaceWellConditions, model, state, ix) where {N, R}
    nstages = length(var.separator_conditions)
    rhoS = vol = missing
    if nstages == 0
        first_time_surface_well_setup!(var, model, state)
    elseif nstages == 1
        rhoS = reference_densities(model.system)
        cond = only(var.separator_conditions)
        rhoS, vol = flash_wellstream_at_surface(model, model.system, state, rhoS, cond)
    else
        rhoS, vol = separator_surface_flash!(var, model, model.system, state)
    end
    x[1] = TopConditions(N, R, density = rhoS, volume_fractions = vol)
end

function first_time_surface_well_setup!(var, model, state)
    sc = physical_representation(model).surface
    push!(var.separator_conditions, sc)
    push!(var.separator_targets, (0, 0))
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

function Base.convert(::Type{TopConditions{N, Float64}}, v::TopConditions{N, <:ForwardDiff.Dual}) where N
    rho = value.(v.density)
    s = value.(v.volume_fractions)
    return TopConditions(N, Float64, density = rho, volume_fractions = s)
end

function Jutul.value(tc::TopConditions{N, <:ForwardDiff.Dual}) where N
    d = value.(tc.density)
    v = value.(tc.volume_fractions)
    return TopConditions(N, Float64, density = d, volume_fractions = v)
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
