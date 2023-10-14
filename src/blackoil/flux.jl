@inline function component_mass_fluxes!(q, face, state, model::StandardBlackOilModel, flux_type, kgrad, upw)
    sys = model.system
    ix = phase_indices(sys)
    (; l, v) = ix
    rhoS = reference_densities(sys)
    kdisc = flux_primitives(face, state, model, flux_type, kgrad, upw)

    # Get the potentials since the flux is more complicated for miscible phases
    potential_difference(phase) = darcy_phase_kgrad_potential(face, phase, state, model, flux_type, kgrad, upw, kdisc)

    b_mob = state.SurfaceVolumeMobilities
    # Water component is the aqueous phase only
    if has_other_phase(sys)
        a = ix.a
        ∇ψ_a = potential_difference(a)
        rhoAS = rhoS[a]
        λb_a = phase_upwind(upw, b_mob, a, ∇ψ_a)
        q_a = rhoAS*λb_a*∇ψ_a
        q = setindex(q, q_a, a)
    end
    # Oil and gas components are always present
    rhoLS = rhoS[l]
    rhoVS = rhoS[v]
    # Oil mobility
    ∇ψ_l = potential_difference(l)
    λb_l = phase_upwind(upw, b_mob, l, ∇ψ_l)
    # Gas mobility
    ∇ψ_v = potential_difference(v)
    λb_v = phase_upwind(upw, b_mob, v, ∇ψ_v)
    if has_vapoil(sys)
        # Rv (vaporized oil) upwinded by vapor potential
        Rv = state.Rv
        f_rv = cell -> @inbounds Rv[cell]
        rv = upwind(upw, f_rv, ∇ψ_v)
        # Final flux = oil phase flux + oil-in-gas flux
        q_l = rhoLS*(λb_l*∇ψ_l + rv*λb_v*∇ψ_v)
    else
        q_l = rhoLS*λb_l*∇ψ_l
    end

    if has_disgas(sys)
        # Rs (solute gas) upwinded by liquid potential
        Rs = state.Rs
        f_rs = cell -> @inbounds Rs[cell]
        rs = upwind(upw, f_rs, ∇ψ_l)
        # Final flux = gas phase flux + gas-in-oil flux
        q_v = rhoVS*(λb_v*∇ψ_v + rs*λb_l*∇ψ_l)
    else
        q_v = rhoVS*λb_v*∇ψ_v
    end

    if haskey(state, :Diffusivities)
        S = state.Saturations
        density = state.PhaseMassDensities
        D = state.Diffusivities
        @inbounds D_l = D[l, face]
        @inbounds D_v = D[v, face]

        if has_disgas(sys)
            X_o = cell -> @inbounds rhoLS/(Rs[cell]*rhoVS)
            X_g = cell -> @inbounds (Rs[cell]*rhoLS)/(Rs[cell]*rhoVS)

            ΔX_o = -gradient(X_o, kgrad)
            ΔX_g = -gradient(X_g, kgrad)

            mass_l = cell -> density[l, cell]*S[l, cell]
            q_l += D_l*upwind(upw, mass_l, ΔX_o)
            q_v += D_l*upwind(upw, mass_l, ΔX_g)
        end

        if has_vapoil(sys)
            Y_o = cell -> @inbounds rhoVS/(Rv[cell]*rhoLS)
            Y_g = cell -> @inbounds (Rv[cell]*rhoVS)/(Rv[cell]*rhoLS)

            ΔY_o = -gradient(Y_o, kgrad)
            ΔY_g = -gradient(Y_g, kgrad)

            mass_v = cell -> density[v, cell]*S[v, cell]
            q_l += D_v*upwind(upw, mass_v, ΔY_o)
            q_v += D_v*upwind(upw, mass_v, ΔY_g)
        end
    end
    q = setindex(q, q_l, l)
    q = setindex(q, q_v, v)
    return q
end

function apply_flow_bc!(acc, q, bc, model::StandardBlackOilModel, state, time)
    mu = state.PhaseViscosities
    b = state.ShrinkageFactors
    kr = state.RelativePermeabilities
    rho = state.PhaseMassDensities
    nph = length(acc)
    @assert size(kr, 1) == nph

    rho_inj = bc.density
    f_inj = bc.fractional_flow
    c = bc.cell
    sys = model.system
    if q > 0
        # Pressure inside is higher than outside, flow out from domain
        phases = phase_indices(sys)
        wat = has_other_phase(sys)
        rhoS = reference_densities(sys)

        if wat
            a, l, v = phases
        else
            l, v = phases
        end

        if wat
            acc[a] += q*rho[a, c]*kr[a, c]/mu[a, c]
        end
        q_l = q_v = 0.0
        q = q*b[l, c]*kr[l, c]/mu[l, c]
        if has_disgas(sys)
            q_v += state.Rs[c]*q
        end
        q_l += q

        q = q*b[v, c]*kr[v, c]/mu[v, c]
        if has_vapoil(sys)
            q_l += state.Rv[c]*q
        end
        q_v += q

        acc[l] += q_l*rhoS[l]
        acc[v] += q_v*rhoS[v]
    else
        # Injection of mass
        λ_t = 0.0
        for ph in eachindex(acc)
            λ_t += kr[ph, c]/mu[ph, c]
        end
        if isnothing(rho_inj)
            # Density not provided, take saturation average from what we have in
            # the inside of the domain
            rho_inj = 0.0
            for ph in 1:nph
                rho_inj += state.Saturations[ph, c]*rho[ph, c]
            end
        end
        if isnothing(f_inj)
            # Fractional flow not provided. We match the mass fraction we
            # observe on the inside.
            total = 0.0
            for ph in 1:nph
                total += state.TotalMasses[ph, c]
            end
            for ph in 1:nph
                F = state.TotalMasses[ph, c]/total
                acc[ph] += q*rho_inj*λ_t*F
            end
        else
            @assert length(f_inj) == nph
            for ph in 1:nph
                F = f_inj[ph]
                acc[ph] += q*rho_inj*λ_t*F
            end
        end
    end
end