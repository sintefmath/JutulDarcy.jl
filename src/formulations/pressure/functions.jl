function store_total_fluxes!(vT, model, state)
    sys = model.system
    nph = number_of_phases(sys)
    N = model.domain.representation.neighborship
    μ = state.PhaseViscosities
    kr = state.RelativePermeabilities
    function mob(ix)
        return 
    end
    for face in eachindex(vT)
        l = N[1, face]
        r = N[2, face]
        tpfa = TPFA(l, r, 1)
        upw = SPU(l, r)
        common = kgrad_common(face, state, model, tpfa)
        v = 0
        for ph in 1:nph
            mob = ix -> kr[ph, ix]/μ[ph, ix]
            q = darcy_phase_kgrad_potential(face, ph, state, model, tpfa, common)
            λ_f = upwind(upw, mob, q)
            v += λ_f * q
        end
        vT[face] = v
    end
end
