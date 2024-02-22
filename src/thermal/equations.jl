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
    q = thermal_heat_flux(face, state, model, grad, upw, flux_type)
    return setindex(Q, q, 1)
end

@inline function Jutul.face_flux!(q_i, face, eq::ConservationLaw{:TotalThermalEnergy, <:Any}, state, model::ThermalModel, dt, flow_disc::PotentialFlow, ldisc)
    # Inner version, for generic flux
    kgrad, upw = ldisc.face_disc(face)
    ft = Jutul.flux_type(eq)
    q = thermal_heat_flux(face, state, model, kgrad, upw, ft)
    return setindex(q_i, q, 1)
end

function thermal_heat_flux(face, state, model, grad, upw, flux_type)
    T = state.Temperature
    H_f = state.FluidEnthalpy
    λ_r = state.RockThermalConductivities
    λ_f = state.FluidThermalConductivities
    S = state.Saturations

    convective_flux = 0
    λ_total = λ_r[face]
    flow_common = kgrad_common(face, state, model, grad)
    for α in 1:number_of_phases(model.system)
        F_α = darcy_phase_mass_flux(face, α, state, model, flux_type, grad, upw, flow_common)
        H_face_α = phase_upwind(upw, H_f, α, F_α)
        convective_flux += H_face_α*F_α
        λ_total += λ_f[α, face]*phase_upwind(upw, S, α, F_α)
    end
    conductive_flux = -λ_total*gradient(T, grad)
    return conductive_flux + convective_flux
end


function Jutul.convergence_criterion(model, storage, eq::ConservationLaw{:TotalThermalEnergy}, eq_s, r; dt = 1.0, update_report = missing)
    a = active_entities(model.domain, Cells())
    E0 = storage.state0.TotalThermalEnergy
    @tullio max e := abs(r[i]) * dt / value(E0[a[i]])
    return (Max = (errors = (e, ), names = ("Energy balance",)), )
end
