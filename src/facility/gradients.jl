function Jutul.vectorization_length(controls_or_limits::Dict, model::FacilityModel, name, variant)
    @assert variant in (:all, :control)
    if variant == :all
        include_temperature = true
        include_mixture_density = true
        include_injection_mixture = true
        include_target = true
    elseif variant == :control
        include_temperature = false
        include_mixture_density = false
        include_injection_mixture = false
        include_target = true
    else
        error("Variant $variant not supported")
    end

    n = 0
    if name == :control
        for (k, v) in pairs(controls_or_limits)
            if v isa DisabledControl
                continue
            else
                if include_target
                    n += 1
                end
                if v isa InjectorControl
                    if include_injection_mixture
                        n += length(v.injection_mixture)
                    end
                    if include_mixture_density
                        n += 1
                    end
                    if include_temperature
                        n += 1
                    end
                else
                    @assert v isa ProducerControl
                end
            end
        end
    else
        error("$name $variant not supported")
    end
    return n
end

function Jutul.vectorize_force!(v, model::FacilityModel, controls_or_limits::Dict, name, variant)
    @assert variant in (:all, :control)
    if variant == :all
        include_temperature = true
        include_mixture_density = true
        include_injection_mixture = true
        include_target = true
    elseif variant == :control
        include_temperature = false
        include_mixture_density = false
        include_injection_mixture = false
        include_target = true
    else
        error("Variant $variant not supported")
    end

    names = Symbol[]
    offset = 0
    if name == :control
        for (wname, ctrl) in pairs(controls_or_limits)
            if ctrl isa DisabledControl
                continue
            else
                if include_target
                    v[offset+1] = ctrl.target.value
                    push!(names, Symbol("target_$wname"))
                    offset += 1
                end
                if ctrl isa InjectorControl
                    if include_injection_mixture
                        for (i, x_i) in enumerate(ctrl.injection_mixture)
                            offset += 1
                            v[offset] = x_i
                            push!(names, Symbol("injection_mixture_$wname$i"))
                        end
                    end
                    if include_mixture_density
                        offset += 1
                        v[offset] = ctrl.mixture_density
                        push!(names, Symbol("mixture_density_$wname"))
                    end
                    if include_temperature
                        offset += 1
                        v[offset] = ctrl.temperature
                        push!(names, Symbol("temperature_$wname"))
                    end
                else
                    @assert ctrl isa ProducerControl
                end
            end
        end
    else
        error("$name not supported")
    end

    return (names = names, )
end
