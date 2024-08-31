function Jutul.select_primary_variables!(vars, ::TransportFormulation, model::TransportModel)
    delete!(vars, :Pressure)
    vars[:TotalSaturation] = TotalSaturation()
end

function Jutul.select_parameters!(vars, ::TransportFormulation, model::TransportModel)
    vars[:Pressure] = Pressure()
    vars[:TotalVolumetricFlux] = TotalVolumetricFlux()
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

@inline function flux_primitives(face, state, model, flux_type::TotalSaturationFlux, tpfa::TPFA, upw)
    trans = state.Transmissibilities
    grav = state.TwoPointGravityDifference
    kr = state.RelativePermeabilities
    mu = state.PhaseViscosities

    @inbounds T_f = trans[face]
    @inbounds gΔz = tpfa.face_sign*grav[face]
    V_t = tpfa.face_sign*state.TotalVolumetricFlux[face]

    ix = phase_indices(model.system)
    l = upw.left
    r = upw.right
    c = upwind_cell(V_t, l, r)
    λ = map(ph -> kr[ph, c]/mu[ph, c], ix)
    λ_t = sum(λ)

    # TODO: Phase upwind and fractional flow here
    return (q = V_t/λ_t, T = T_f, gdz = gΔz, V_t = V_t, λ = λ)
end

@inline function darcy_phase_kgrad_potential(face, phase, state, model, flux_type::TotalSaturationFlux, tpfa::TPFA{T}, upw, common = flux_primitives(face, state, model, flux_type, upw, tpfa)) where T
    ρ = state.PhaseMassDensities
    pc, ref_index = capillary_pressure(model, state)
    V_t, T_f, gΔz = common

    l = tpfa.left
    r = tpfa.right
    pc::Nothing
    @assert gΔz == 0
    # Δpc = capillary_gradient(pc, l, r, phase, ref_index)
    # @inbounds ρ_c = ρ[phase, l]
    # @inbounds ρ_i = ρ[phase, r]
    ## ρ_avg = 0.5*(ρ_i + ρ_c)
    # q = -T_f*(∇p + Δpc + gΔz*ρ_avg)
    # error()
    q = V_t
    return q
end
