@inline function Jutul.face_flux!(Q, left, right, face, face_sign, eq::ConservationLaw{:TracerMasses, <:Any}, state, model, dt, flow_disc::TwoPointPotentialFlowHardCoded)
    kgrad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    flux_type = Jutul.flux_type(eq)
    return tracer_flux(Q, face, state, model, kgrad, upw, flux_type)
end

@inline function Jutul.face_flux!(q_i, face, eq::ConservationLaw{:TracerMasses, <:Any}, state, model, dt, flow_disc::PotentialFlow, ldisc)
    kgrad, upw = ldisc.face_disc(face)
    ft = Jutul.flux_type(eq)
    return tracer_flux(Q, face, state, model, kgrad, upw, ft)
end

function tracer_flux(Q, face, state, model, kgrad, upw, ft::TracerFluxType)
    phases = tuple(1:number_of_phases(model.system)...)
    flow_common = kgrad_common(face, state, model, kgrad)
    phase_mass_fluxes = map(α -> darcy_phase_mass_flux(face, α, state, model, ft, kgrad, upw, flow_common), phases)

    C = state.TracerConcentrations
    tracers = ft.tracers
    N = length(tracers)
    T = eltype(Q)
    for i in 1:N
        tracer = tracers[i]
        v = zero(T)
        for phase in tracer_phase_indices(tracer)
            q_ph = phase_mass_fluxes[phase]
            C_iface = phase_upwind(upw, C, phase, q_ph)
            v += C_iface*q_ph
        end
        Q = setindex(Q, v, i)
    end
    return Q
end

function Jutul.convergence_criterion(model, storage, eq::ConservationLaw{:TracerMasses}, eq_s, r; dt = 1.0, update_report = missing)
    a = active_entities(model.domain, Cells())
    V = storage.state0.FluidVolume
    nt = size(r, 1)
    e = zeros(nt)
    for c in a
        scale = dt/V[c]
        for t in 1:nt
            v = abs(r[t, c])*scale
            e[t] = max(e[t], v)
        end
    end
    cnames = model.equations[:tracers].flux_type.names
    return (Max = (errors = e, names = cnames), )
end
