"""
    component_mass_fluxes!(q, face, state, model::StandardBlackOilModel, flux_type, kgrad, upw)

Compute the component mass fluxes for a given face in a black oil model.

# Arguments
- `q`: Output array to store the computed mass fluxes.
- `face`: The face for which the fluxes are being computed.
- `state`: The current state of the system.
- `model::StandardBlackOilModel`: The black oil model being used.
- `flux_type`: The type of flux calculation to be performed.
- `kgrad`: Gradient operator for K grad p.
- `upw`: Upwind scheme to be used for the flux calculation.

# Returns
- This function modifies the `q` array in place with the computed mass fluxes.
"""
@inline function component_mass_fluxes!(q, face, state, model::StandardBlackOilModel, flux_type, kgrad, upw)
    sys = model.system
    ix = phase_indices(sys)
    (; l, v) = ix
    rhoS = reference_densities(sys)
    # Get the potentials since the flux is more complicated for miscible phases
    potential_differences = darcy_permeability_potential_differences(face, state, model, flux_type, kgrad, upw)

    b_mob = state.SurfaceVolumeMobilities
    # Water component is the aqueous phase only
    if has_other_phase(sys)
        a = ix.a
        ∇ψ_a = potential_differences[a]
        rhoAS = rhoS[a]
        λb_a = phase_upwind(upw, b_mob, a, ∇ψ_a)
        q_a = rhoAS*λb_a*∇ψ_a
        q = setindex(q, q_a, a)
    end
    # Oil and gas components are always present
    rhoLS = rhoS[l]
    rhoVS = rhoS[v]
    # Oil mobility
    ∇ψ_l = potential_differences[l]
    λb_l = phase_upwind(upw, b_mob, l, ∇ψ_l)
    # Gas mobility
    ∇ψ_v = potential_differences[v]
    λb_v = phase_upwind(upw, b_mob, v, ∇ψ_v)
    if has_vapoil(sys)
        # Rv (vaporized oil) upwinded by vapor potential
        Rv = state.Rv
        rv = upwind(upw, Rv, ∇ψ_v)
        # Final flux = oil phase flux + oil-in-gas flux
        q_l = rhoLS*(λb_l*∇ψ_l + rv*λb_v*∇ψ_v)
    else
        q_l = rhoLS*λb_l*∇ψ_l
    end

    if has_disgas(sys)
        # Rs (solute gas) upwinded by liquid potential
        Rs = state.Rs
        rs = upwind(upw, Rs, ∇ψ_l)
        # Final flux = gas phase flux + gas-in-oil flux
        q_v = rhoVS*(λb_v*∇ψ_v + rs*λb_l*∇ψ_l)
    else
        q_v = rhoVS*λb_v*∇ψ_v
    end

    if haskey(state, :Diffusivities)
        S = state.Saturations
        density = state.PhaseMassDensities
        D = state.Diffusivities
        if has_disgas(sys)
            qo_diffusive_l, qo_diffusive_v = blackoil_diffusion(Rs, S, density, rhoLS, rhoVS, face, D, l, kgrad, upw)
            q_l += qo_diffusive_l
            q_v += qo_diffusive_v
        end

        if has_vapoil(sys)
            qg_diffusive_v, qg_diffusive_l = blackoil_diffusion(Rv, S, density, rhoVS, rhoLS, face, D, v, kgrad, upw)
            q_l += qg_diffusive_l
            q_v += qg_diffusive_v
        end
    end
    q = setindex(q, q_l, l)
    q = setindex(q, q_v, v)
    return q
end

function blackoil_diffusion(R, S, density, rhoS_self, rhoS_dissolved, face, D, α, kgrad, upw)
    @inbounds D_α = D[α, face]
    X_self = cell -> black_oil_phase_mass_fraction(rhoS_self, rhoS_dissolved, R, cell)
    # Two components: 1 - X_l - (1 - X_r) = - X_l + X_r = -(X_l - X_r) = ΔX
    ΔX_self = -gradient(X_self, kgrad)
    ΔX_other = -ΔX_self

    T = typeof(ΔX_self)
    mass_l = cell -> density[α, cell]*S[α, cell]
    # TODO: Upwind or average here? Maybe doesn't matter, should be in
    # parabolic limit for diffusion
    # q_l += D_l*upwind(upw, mass_l, ΔX_o)*ΔX_o
    # q_v += D_l*upwind(upw, mass_l, ΔX_g)*ΔX_g

    diffused_mass = D_α*face_average(mass_l, kgrad)
    diff_self = convert(T, diffused_mass*ΔX_self)
    diff_dissolved = convert(T, diffused_mass*ΔX_other)
    return (diff_self::T, diff_dissolved::T)::Tuple{T, T}
end

@inline function black_oil_phase_mass_fraction(rhoLS, rhoVS, Rs, cell)
    # TODO: Should have molar weights here maybe, but not part of standard input
    @inbounds rs = Rs[cell]
    if rs < 1e-10
        v = one(rs)
    else
        v = rhoLS/(rhoLS + rs*rhoVS)
    end
    return v
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
