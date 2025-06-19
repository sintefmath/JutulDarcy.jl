@inline function Jutul.face_flux!(Q, left, right, face, face_sign, eq::ConservationLaw{:TracerMasses, <:Any}, state, model, dt, flow_disc::TwoPointPotentialFlowHardCoded)
    kgrad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    flux_type = Jutul.flux_type(eq)
    return tracer_flux(Q, face, state, model, kgrad, upw, flux_type)
end

@inline function Jutul.face_flux!(Q, face, eq::ConservationLaw{:TracerMasses, <:Any}, state, model, dt, flow_disc::PotentialFlow, ldisc)
    kgrad, upw = ldisc.face_disc(face)
    ft = Jutul.flux_type(eq)
    return tracer_flux(Q, face, state, model, kgrad, upw, ft)
end

function tracer_flux(Q, face, state, model, kgrad, upw, ft::TracerFluxType)
    phase_mass_fluxes = JutulDarcy.darcy_phase_mass_fluxes(face, state, model, ft, kgrad, upw)

    tracers = ft.tracers
    N = length(tracers)
    T = eltype(Q)
    for i in 1:N
        tracer = tracers[i]
        v = tracer_total_mass_flux(tracer, model, state, phase_mass_fluxes, i, upw)
        Q = setindex(Q, v, i)
    end
    return Q
end

function tracer_scale(model, tracer::AbstractTracer)
    # Assume density of roughly 100.0 This function reduces the convergence
    # criterion by this factor. Higher values means slacker criterion.
    return 100.0
end

function Jutul.convergence_criterion(model, storage, eq::ConservationLaw{:TracerMasses}, eq_s, r; dt = 1.0, update_report = missing)
    a = active_entities(model.domain, Cells())
    V = storage.state0.FluidVolume
    nt = size(r, 1)
    e = zeros(nt)
    tracers = eq.flux_type.tracers
    tscale = map(t -> tracer_scale(model, t), tracers)
    for c in a
        scale = dt/V[c]
        for t in 1:nt
            v = abs(r[t, c])*scale/tscale[t]
            e[t] = max(e[t], v)
        end
    end
    cnames = model.equations[:tracers].flux_type.names
    return (Max = (errors = e, names = cnames), )
end
