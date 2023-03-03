@inline function component_mass_fluxes!(q, face, state, model::StandardBlackOilModel, kgrad, upw)
    sys = model.system
    a, l, v = phase_indices(sys)
    rhoAS, rhoLS, rhoVS = reference_densities(sys)
    kdisc = kgrad_common(face, state, model, kgrad)

    # Get the potentials since the flux is more complicated for miscible phases
    potential(phase) = darcy_phase_kgrad_potential(face, phase, state, model, kgrad, kdisc)

    b_mob = state.SurfaceVolumeMobilities
    # Water component is the aqueous phase only
    ψ_a = potential(a)
    λb_a = phase_upwind(upw, b_mob, a, ψ_a)
    q_a = rhoAS*λb_a*ψ_a

    # Oil mobility
    ψ_l = potential(l)
    λb_l = phase_upwind(upw, b_mob, l, ψ_l)
    # Gas mobility
    ψ_v = potential(v)
    λb_v = phase_upwind(upw, b_mob, v, ψ_v)

    # Rv (vaporized oil) upwinded by vapor potential
    Rv = state.Rv
    f_rv = cell -> @inbounds Rv[cell]
    rv = upwind(upw, f_rv, ψ_v)
    # Final flux = oil phase flux + oil-in-gas flux
    q_l = rhoLS*(λb_l*ψ_l + rv*λb_v*ψ_v)

    # Rs (solute gas) upwinded by liquid potential
    Rs = state.Rs
    f_rs = cell -> @inbounds Rs[cell]
    rs = upwind(upw, f_rs, ψ_l)
    # Final flux = gas phase flux + gas-in-oil flux
    q_v = rhoVS*(λb_v*ψ_v + rs*λb_l*ψ_l)

    q = setindex(q, q_a, a)
    q = setindex(q, q_l, l)
    q = setindex(q, q_v, v)
    return q
end

@inline function component_mass_fluxes!(q, face, state, model::VapoilBlackOilModel, kgrad, upw)
    sys = model.system
    a, l, v = phase_indices(sys)
    rhoAS, rhoLS, rhoVS = reference_densities(sys)
    kdisc = kgrad_common(face, state, model, kgrad)

    # Get the potentials since the flux is more complicated for miscible phases
    potential(phase) = darcy_phase_kgrad_potential(face, phase, state, model, kgrad, kdisc)

    b_mob = state.SurfaceVolumeMobilities
    # Water component is the aqueous phase only
    ψ_a = potential(a)
    λb_a = phase_upwind(upw, b_mob, a, ψ_a)
    q_a = rhoAS*λb_a*ψ_a

    # Oil mobility
    ψ_l = potential(l)
    λb_l = phase_upwind(upw, b_mob, l, ψ_l)
    # Gas mobility
    ψ_v = potential(v)
    λb_v = phase_upwind(upw, b_mob, v, ψ_v)

    # Rv (vaporized oil) upwinded by vapor potential
    Rv = state.Rv
    f_rv = cell -> @inbounds Rv[cell]
    rv = upwind(upw, f_rv, ψ_v)
    # Final flux = oil phase flux + oil-in-gas flux
    q_l = rhoLS*(λb_l*ψ_l + rv*λb_v*ψ_v)
    # Vapor flux is simple
    q_v = rhoVS*λb_v*ψ_v

    q = setindex(q, q_a, a)
    q = setindex(q, q_l, l)
    q = setindex(q, q_v, v)
    return q
end

@inline function component_mass_fluxes!(q, face, state, model::DisgasBlackOilModel, kgrad, upw)
    sys = model.system
    has_water = has_other_phase(sys)
    if has_water
        a, l, v = phase_indices(sys)
        rhoAS, rhoLS, rhoVS = reference_densities(sys)
    else
        l, v = phase_indices(sys)
        rhoLS, rhoVS = reference_densities(sys)
    end
    kdisc = kgrad_common(face, state, model, kgrad)

    # Get the potentials since the flux is more complicated for miscible phases
    potential(phase) = darcy_phase_kgrad_potential(face, phase, state, model, kgrad, kdisc)

    b_mob = state.SurfaceVolumeMobilities
    if has_water
        # Water component is the aqueous phase only
        ψ_a = potential(a)
        λb_a = phase_upwind(upw, b_mob, a, ψ_a)
        q_a = rhoAS*λb_a*ψ_a
        q = setindex(q, q_a, a)
    end
    # Oil component is the oil phase only
    ψ_l = potential(l)
    λb_l = phase_upwind(upw, b_mob, l, ψ_l)
    q_l = rhoLS*λb_l*ψ_l
    # Gas mobility
    ψ_v = potential(v)
    λb_v = phase_upwind(upw, b_mob, v, ψ_v)
    # Rs (solute gas) upwinded by liquid potential
    Rs = state.Rs
    f_rs = cell -> @inbounds Rs[cell]
    rs = upwind(upw, f_rs, ψ_l)
    # Final flux = gas phase flux + gas-in-oil flux
    q_v = rhoVS*(λb_v*ψ_v + rs*λb_l*ψ_l)

    q = setindex(q, q_l, l)
    q = setindex(q, q_v, v)
    return q
end
