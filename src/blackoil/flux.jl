@inline function component_mass_fluxes!(q, face, state, model::SimulationModel{<:Any, <:BlackOilSystem, <:Any, <:Any}, kgrad, upw)
    sys = model.system
    a, l, v = phase_indices(sys)
    rhoS = reference_densities(sys)
    # Aqueous phase is straightforward, just the mass flux for phase
    q_a = darcy_phase_mass_flux(face, a, state, model, kgrad, upw)
    # Get the potentials since the flux is more complicated for miscible phases
    ψ_l = darcy_phase_kgrad_potential(face, l, state, model, kgrad)
    ψ_v = darcy_phase_kgrad_potential(face, v, state, model, kgrad)

    kr = state.RelativePermeabilities
    μ = state.PhaseViscosities
    b = state.ShrinkageFactors
    Rs = state.Rs

    # Oil
    f_bl = cell -> b_mob(b, kr, μ, l, cell)
    λb_l = upwind(upw, f_bl, ψ_l)
    q_l = rhoS[l]*λb_l*ψ_l
    # Gas mobility
    f_bv = cell -> b_mob(b, kr, μ, v, cell)
    λb_v = upwind(upw, f_bv, ψ_v)
    # Rs (solute gas) upwinded by liquid potential
    f_rs = cell -> @inbounds Rs[cell]
    rs = upwind(upw, f_rs, ψ_l)
    # Final flux = gas phase flux + gas-in-oil flux
    q_v = (λb_v*ψ_v + rs*λb_l*ψ_l)*rhoS[v]

    q = setindex(q, q_a, a)
    q = setindex(q, q_l, l)
    q = setindex(q, q_v, v)
    return q
end

function b_mob(b, kr, μ, ph, c)
    @inbounds λ = kr[ph, c]/μ[ph, c]
    @inbounds b_f = b[ph, c]
    return λ*b_f
end
