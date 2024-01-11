function select_equations!(eqs, sys::ThermalSystem, model::SimulationModel)
    disc = model.domain.discretizations.heat_flow
    eqs[:energy_conservation] = ConservationLaw(disc, :TotalThermalEnergy, 1)
end

@inline function Jutul.face_flux!(Q, left, right, face, face_sign, eq::ConservationLaw{:TotalThermalEnergy, <:Any}, state, model::ThermalModel, dt, flow_disc::TwoPointPotentialFlowHardCoded)
    # Specific version for tpfa flux
    # TODO: Add general version for thermal
    grad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    # TODO: This could be inconsistent if we use flux_type for something.
    flux_type = Jutul.flux_type(eq)
    q = thermal_heat_flux!(face, state, model, grad, upw, flux_type)
    return setindex(Q, q, 1)
end

function thermal_heat_flux!(face, state, model, grad, upw, flux_type)
    T = state.Temperature
    H_f = state.FluidEnthalpy
    λ_r = state.RockThermalConductivities[face]
    λ_f = state.FluidThermalConductivities[face]

    conductive_flux = -(λ_r + λ_f)*gradient(T, grad)
    convective_flux = zero(conductive_flux)

    flow_common = kgrad_common(face, state, model, grad)
    for α in 1:number_of_phases(model.system)
        F_α = darcy_phase_mass_flux(face, α, state, model, flux_type, grad, upw, flow_common)
        H_face_α = phase_upwind(upw, H_f, α, F_α)
        convective_flux += H_face_α*F_α
    end
    return conductive_flux + convective_flux
end
