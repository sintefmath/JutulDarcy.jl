function store_total_fluxes!(vT, model, state)
    sys = model.system
    nph = number_of_phases(sys)
    N = model.domain.representation.neighborship
    @assert length(vT) == size(N, 2) "vT must have the same length (was $(length(vT))) as the number of faces $(size(N, 2))"
    μ = state.PhaseViscosities
    kr = state.RelativePermeabilities
    for face in eachindex(vT)
        l = N[1, face]
        r = N[2, face]
        tpfa = TPFA(l, r, 1)
        upw = SPU(l, r)
        f_t = Jutul.DefaultFlux()
        # TODO: This assumes the default discretizations.
        v_phases = darcy_phase_volume_fluxes(face, state, model, f_t, tpfa, upw)
        vT[face] = sum(v_phases)
    end
    return vT
end

function store_total_fluxes(model, state::AbstractDict)
    return store_total_fluxes(model, JutulStorage(state))
end

function store_total_fluxes(model, state::Union{JutulStorage, NamedTuple})
    nf = number_of_faces(model.data_domain)
    vT = zeros(nf)
    return store_total_fluxes!(vT, model, state)
end

function store_phase_fluxes!(v_phases, model, state)
    sys = model.system
    nph = number_of_phases(sys)
    N = model.domain.representation.neighborship

    @assert size(v_phases, 1) == nph
    @assert size(v_phases, 2) == size(N, 2)
    for face in axes(v_phases, 2)
        l = N[1, face]
        r = N[2, face]
        tpfa = TPFA(l, r, 1)
        upw = SPU(l, r)
        f_t = Jutul.DefaultFlux()
        common = flux_primitives(face, state, model, f_t, tpfa, upw)
        for ph in 1:nph
            mob = ix -> kr[ph, ix]/μ[ph, ix]
            q = darcy_phase_kgrad_potential(face, ph, state, model, f_t, tpfa, upw, common)
            λ_f = upwind(upw, mob, q)
            v_phases[ph, face] = λ_f * q
        end
    end
    return v_phases
end

function store_phase_fluxes(model, state)
    nf = number_of_faces(model.data_domain)
    nph = number_of_phases(model.system)
    v = zeros(nph, nf)
    return store_phase_fluxes!(v, model, state)
end
