struct TracerConcentrations{T} <: Jutul.VectorVariables
    flux_type::T
end

Jutul.values_per_entity(model, t::TracerConcentrations) = number_of_tracers(t.flux_type)
Jutul.minimum_value(t::TracerConcentrations) = 0.0

struct TracerMasses{T} <: Jutul.VectorVariables
    flux_type::T
end

Jutul.values_per_entity(model, t::TracerMasses) = number_of_tracers(t.flux_type)

Jutul.@jutul_secondary function update_tracer_masses!(tracer_mass, tm::TracerMasses, model, TracerConcentrations, Saturations, PhaseMassDensities, FluidVolume, ix)
    tracers = tm.flux_type.tracers
    N = length(tracers)
    T = eltype(tracer_mass)
    for cell in ix
        vol = FluidVolume[cell]
        for i in 1:N
            tracer = tracers[i]
            resident_mass = zero(T)
            c = TracerConcentrations[i, cell]
            for phase in tracer_phase_indices(tracer)
                S = Saturations[phase, cell]
                den = PhaseMassDensities[phase, cell]
                resident_mass += S*den
            end
            if resident_mass <= TRACER_TOL
                v = Jutul.replace_value(c, 0.0)
                # v = c*value(resident_mass)
            else
                v = c*resident_mass
            end
            tracer_mass[i, cell] = v*vol
        end
    end
    return tracer_mass
end
