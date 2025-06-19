function Jutul.select_primary_variables!(vars, ::TransportFormulation, model::TransportModel)
    delete!(vars, :Pressure)
    vars[:TotalSaturation] = TotalSaturation()
end

function Jutul.select_parameters!(vars, ::TransportFormulation, model::TransportModel)
    vars[:Pressure] = Pressure()
    if JutulDarcy.model_or_domain_is_well(model)
        vars[:PerforationTotalVolumetricFlux] = PerforationTotalVolumetricFlux()
    else
        vars[:TotalVolumetricFlux] = TotalVolumetricFlux()
    end
end

function Jutul.select_equations!(eqs, ::TransportFormulation, model::TransportModel)
    lbl = :mass_conservation
    @assert haskey(eqs, lbl)
    eq = eqs[lbl]

    s = conserved_symbol(eq)
    @assert s == :TotalMasses
    N = number_of_equations_per_entity(model, eq)
    eqs[lbl] = ConservationLaw(eq.flow_discretization, s, N, flux = TotalSaturationFlux())
end

@inline function JutulDarcy.darcy_permeability_potential_differences(
        face,
        state,
        model,
        flux_type::TotalSaturationFlux,
        kgrad::TPFA,
        upw,
        phases = JutulDarcy.eachphase(model.system)
    )
    V_t = kgrad.face_sign*state.TotalVolumetricFlux[face]
    N = length(phases)
    @assert N == number_of_phases(model.system)
    l, r = Jutul.cell_pair(upw)

    left_mob = map(phase -> state.PhaseMobilities[phase, l], phases)
    right_mob = map(phase -> state.PhaseMobilities[phase, r], phases)

    T_f = effective_transmissibility(state, face, kgrad)
    gΔz = effective_gravity_difference(state, face, kgrad)
    pc, ref_index = capillary_pressure(model, state)

    @inline function bouyancy_and_capillary_term(phase)
        Δpc = capillary_gradient(pc, kgrad, phase, ref_index)
        ρ_avg = face_average_density(model, state, kgrad, phase)
        return -(gΔz*ρ_avg + Δpc)
        # return -T_f*(∇p + Δpc + gΔz*ρ_avg)
    end

    G = map(bouyancy_and_capillary_term, phases)

    # ∇p = pressure_gradient(state, kgrad)

    # @inline function phase_pot(phase)
    #     Δpc = capillary_gradient(pc, kgrad, phase, ref_index)
    #     ρ_avg = face_average_density(model, state, kgrad, phase)
    #     return -T_f*(∇p + Δpc + gΔz*ρ_avg)
    # end
    if false
        mob_1 = JutulDarcy.phase_upwind(upw, state.PhaseMobilities, 1, V_t)
        mob_2 = JutulDarcy.phase_upwind(upw, state.PhaseMobilities, 2, V_t)
        mob_t = mob_1 + mob_2
        @info "???" V_t mob_1 mob_2
        return (1/mob_t*V_t, 1/mob_t*V_t)
    end
    return phase_potential_upwind_potential_differences(V_t, T_f, G, left_mob, right_mob)
end


# @inline function flux_primitives(face, state, model, flux_type::TotalSaturationFlux, tpfa::TPFA, upw)
#     trans = state.Transmissibilities
#     grav = state.TwoPointGravityDifference
#     kr = state.RelativePermeabilities
#     mu = state.PhaseViscosities

#     @inbounds T_f = trans[face]
#     @inbounds gΔz = tpfa.face_sign*grav[face]
#     V_t = tpfa.face_sign*state.TotalVolumetricFlux[face]

#     ix = phase_indices(model.system)
#     l = upw.left
#     r = upw.right
#     c = upwind_cell(V_t, l, r)
#     λ = map(ph -> kr[ph, c]/mu[ph, c], ix)
#     λ_t = sum(λ)

#     # TODO: Phase upwind and fractional flow here
#     return (q = V_t/λ_t, T = T_f, gdz = gΔz, V_t = V_t, λ = λ)
# end

# @inline function darcy_phase_kgrad_potential(face, phase, state, model, flux_type::TotalSaturationFlux, tpfa::TPFA{T}, upw, common = flux_primitives(face, state, model, flux_type, upw, tpfa)) where T
#     ρ = state.PhaseMassDensities
#     pc, ref_index = capillary_pressure(model, state)
#     V_t, T_f, gΔz = common

#     l = tpfa.left
#     r = tpfa.right
#     pc::Nothing
#     @assert gΔz == 0
#     # Δpc = capillary_gradient(pc, l, r, phase, ref_index)
#     # @inbounds ρ_c = ρ[phase, l]
#     # @inbounds ρ_i = ρ[phase, r]
#     ## ρ_avg = 0.5*(ρ_i + ρ_c)
#     # q = -T_f*(∇p + Δpc + gΔz*ρ_avg)
#     # error()
#     q = V_t
#     return q
# end

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t,
        state_s, state0_s, 
        model_t::TransportModel, model_s::TransportModel,
        ct::ReservoirFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i)
    )
    error("This function is not implemented for TransportModel.")
    sys = model_t.system
    rhoS = reference_densities(sys)
    conn = cross_term_perforation_get_conn(ct, i, state_s, state_t)
    # Call smaller interface that is easy to specialize
    if haskey(state_s, :MassFractions)
        @inbounds simple_well_perforation_flux!(out, sys, state_t, state_s, rhoS, conn)
    else
        @inbounds multisegment_well_perforation_flux!(out, sys, state_t, state_s, rhoS, conn)
    end
end
