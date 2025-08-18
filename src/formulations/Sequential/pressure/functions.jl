function store_phase_fluxes(model, state, is_mass = false)
    nf = number_of_faces(model.data_domain)
    nph = number_of_phases(model.system)
    v = zeros(nph, nf)
    return store_phase_fluxes!(v, model, state, is_mass)
end

function store_total_fluxes(model, state::AbstractDict, is_mass = false)
    return store_total_fluxes(model, JutulStorage(state), is_mass)
end

function store_total_fluxes!(vT, model, state, is_mass::Bool = false)
    N = model.domain.representation.neighborship
    @assert length(vT) == size(N, 2) "vT must have the same length (was $(length(vT))) as the number of faces $(size(N, 2))"
    for face in eachindex(vT)
        v_phases = store_flux_helper(face, N, state, model, is_mass)
        vT[face] = sum(v_phases)
    end
    return vT
end

function store_total_fluxes(model, state::Union{JutulStorage, NamedTuple}, is_mass = false)
    nf = number_of_faces(model.data_domain)
    vT = zeros(nf)
    return store_total_fluxes!(vT, model, state, is_mass)
end

function store_phase_fluxes!(v_phases, model, state, is_mass::Bool = false)
    sys = model.system
    nph = number_of_phases(sys)
    N = model.domain.representation.neighborship

    @assert size(v_phases, 1) == nph
    @assert size(v_phases, 2) == size(N, 2)
    for face in axes(v_phases, 2)
        v_face = store_flux_helper(face, N, state, model, is_mass)
        for ph in 1:nph
            v_phases[ph, face] = v_face[ph]
        end
    end
    return v_phases
end

function store_flux_helper(face, N, state, model, is_mass::Bool)
    l = N[1, face]
    r = N[2, face]
    # TODO: This assumes the default discretizations.
    # This should be generalized to allow for different flux types.
    tpfa = TPFA(l, r, 1)
    upw = SPU(l, r)
    f_t = Jutul.DefaultFlux()

    if is_mass
        v_face = JutulDarcy.darcy_phase_mass_fluxes(face, state, model, f_t, tpfa, upw)
    else
        v_face = JutulDarcy.darcy_phase_volume_fluxes(face, state, model, f_t, tpfa, upw)
    end
    return v_face
end