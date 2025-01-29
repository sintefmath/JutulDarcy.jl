struct TracerConcentrations{T} <: Jutul.VectorVariables
    flux_type::T
end

Jutul.values_per_entity(model, t::TracerConcentrations) = number_of_tracers(t.flux_type)
Jutul.minimum_value(t::TracerConcentrations) = 0.0

struct TracerMasses{T} <: Jutul.VectorVariables
    flux_type::T
end

Jutul.values_per_entity(model, t::TracerMasses) = number_of_tracers(t.flux_type)

function Jutul.update_secondary_variable!(tracer_mass, tm::TracerMasses, model, state, ix)
    TracerConcentrations = state.TracerConcentrations
    Saturations = state.Saturations
    PhaseMassDensities = state.PhaseMassDensities
    FluidVolume = state.FluidVolume

    tracers = tm.flux_type.tracers
    N = length(tracers)
    T = eltype(tracer_mass)
    for cell in ix
        vol = FluidVolume[cell]
        for i in 1:N
            tracer = tracers[i]
            resident_mass_density = zero(T)
            c = TracerConcentrations[i, cell]
            for phase in tracer_phase_indices(tracer)
                S = Saturations[phase, cell]
                den = PhaseMassDensities[phase, cell]
                resident_mass_density += S*den
            end
            tracer_mass[i, cell] = tracer_total_mass_outer(tracer, model, state, c, resident_mass_density, vol, cell, i)
        end
    end
    return tracer_mass
end

function Jutul.get_dependencies(var::TracerMasses, model)
    out = Symbol[
        :TracerConcentrations,
        :Saturations,
        :PhaseMassDensities,
        :FluidVolume
    ]
    for tracer in var.flux_type.tracers
        for dep in Jutul.get_dependencies(tracer, model)
            push!(out, dep)
        end
    end
    return unique!(out)
end
