function facility_vectorization(variant::Symbol)
    @assert variant in (:all, :control)
    if variant == :all
        # Default behavior
        include_temperature = true
        include_mixture_density = false
        include_injection_mixture = false
        include_target = true
        include_limits = true
    elseif variant == :control
        include_temperature = false
        include_mixture_density = false
        include_injection_mixture = false
        include_target = true
        include_limits = false
    elseif variant == :extended
        include_temperature = true
        include_mixture_density = true
        include_injection_mixture = true
        include_target = true
        include_limits = true
    else
        error("Variant $variant not supported")
    end
    return (
            temperature = include_temperature,
            mixture_density = include_mixture_density,
            injection_mixture = include_injection_mixture,
            target = include_target,
            limits = include_limits
        )
end

function Jutul.vectorization_length(controls_or_limits::AbstractDict, model::FacilityModel, name, variant)
    supp = facility_vectorization(variant)
    n = 0
    if name == :control
        for (k, v) in pairs(controls_or_limits)
            inner_v = v isa GroupControl ? v.well_control : v
            if inner_v isa DisabledControl
                continue
            else
                if supp.target
                    n += 1
                end
                if inner_v isa InjectorControl
                    if supp.injection_mixture
                        n += length(inner_v.injection_mixture)
                    end
                    if supp.mixture_density
                        n += 1
                    end
                    if supp.temperature
                        n += 1
                    end
                else
                    @assert inner_v isa ProducerControl
                end
            end
        end
    elseif name == :limits
        if supp.limits
            for (k, v) in pairs(controls_or_limits)
                if !isnothing(v)
                    for (lim_k, lim_v) in pairs(v)
                        n += 1
                    end
                end
            end
        end
    else
        error("$name $variant not supported")
    end
    return n
end

function Jutul.vectorize_force!(v, model::FacilityModel, controls_or_limits::AbstractDict, name, variant)
    supp = facility_vectorization(variant)
    names = Symbol[]
    offset = 0
    if name == :control
        for (wname, ctrl) in pairs(controls_or_limits)
            inner_ctrl = ctrl isa GroupControl ? ctrl.well_control : ctrl
            if inner_ctrl isa DisabledControl
                continue
            else
                if supp.target
                    v[offset+1] = inner_ctrl.target.value
                    push!(names, Symbol("target_$wname"))
                    offset += 1
                end
                if inner_ctrl isa InjectorControl
                    if supp.injection_mixture
                        for (i, x_i) in enumerate(inner_ctrl.injection_mixture)
                            offset += 1
                            v[offset] = x_i
                            push!(names, Symbol("injection_mixture_$(wname)_$i"))
                        end
                    end
                    if supp.mixture_density
                        offset += 1
                        v[offset] = inner_ctrl.mixture_density
                        push!(names, Symbol("mixture_density_$wname"))
                    end
                    if supp.temperature
                        offset += 1
                        v[offset] = inner_ctrl.temperature
                        push!(names, Symbol("temperature_$wname"))
                    end
                else
                    @assert inner_ctrl isa ProducerControl
                end
            end
        end
    elseif name == :limits
        if supp.limits
            for (k, limdict) in pairs(controls_or_limits)
                if isnothing(limdict)
                    continue
                end
                for (lim_k, lim_v) in pairs(limdict)
                    offset += 1
                    v[offset] = lim_v
                    push!(names, Symbol("limit_$k$lim_k"))
                end
            end
        end
    else
        error("$name $variant not supported")
    end

    return (names = names, )
end

function Jutul.devectorize_force(control_or_limits::Tcl, model::FacilityModel, X, meta, name, variant) where Tcl
    control_or_limits::AbstractDict
    supp = facility_vectorization(variant)
    offset = 0
    out = Tcl()
    T = eltype(X)
    if name == :control
        for (wname, ctrl) in pairs(control_or_limits)
            inner_ctrl = ctrl isa GroupControl ? ctrl.well_control : ctrl
            if inner_ctrl isa DisabledControl
                out[wname] = ctrl
            else
                if supp.target
                    val = X[offset+1]
                    if inner_ctrl.target isa ReinjectionTarget
                        target = deepcopy(inner_ctrl.target)
                        target.value = val
                    else
                        Tt = Base.typename(typeof(inner_ctrl.target)).wrapper
                        # TODO: Ugly hack...
                        target = Tt(val)
                    end
                    offset += 1
                else
                    target = T(inner_ctrl.target)
                end
                if inner_ctrl isa InjectorControl
                    if supp.injection_mixture
                        nm = length(inner_ctrl.injection_mixture)
                        mixture = X[offset+1:offset+nm]
                        offset += nm
                    else
                        mixture = T.(inner_ctrl.injection_mixture)
                    end
                    if supp.mixture_density
                        density = X[offset+1]
                        offset += 1
                    else
                        density = T(inner_ctrl.mixture_density)
                    end
                    if supp.temperature
                        temp = X[offset+1]
                        offset += 1
                    else
                        temp = T(inner_ctrl.temperature)
                    end
                    new_inner = InjectorControl(target, mixture,
                        density = density,
                        temperature = temp,
                        factor = T(inner_ctrl.factor),
                        enthalpy = inner_ctrl.enthalpy,
                        tracers = inner_ctrl.tracers,
                        phases = inner_ctrl.phases,
                        check = false
                    )
                    if ctrl isa GroupControl
                        out[wname] = GroupControl(new_inner, ctrl.group, allocation_factor = ctrl.allocation_factor)
                    else
                        out[wname] = new_inner
                    end
                else
                    @assert inner_ctrl isa ProducerControl
                    new_inner = ProducerControl(target, check = false)
                    if ctrl isa GroupControl
                        out[wname] = GroupControl(new_inner, ctrl.group, allocation_factor = ctrl.allocation_factor)
                    else
                        out[wname] = new_inner
                    end
                end
            end
        end
    elseif name == :limits
        if supp.limits
            for (k, limdict) in pairs(control_or_limits)
                new_limdict = OrderedDict()
                if isnothing(limdict)
                    continue
                end
                for (lim_k, lim_v) in pairs(limdict)
                    offset += 1
                    new_limdict[lim_k] = X[offset]
                end
                out[k] = (; pairs(new_limdict)...)
            end
        end
    else
        error("$name $variant not supported")
    end
    return out
end