const EnergyConservationLaw = Union{
    ConservationLaw{:TotalThermalEnergy, <:Any},
    ConservationLaw{:TotalEnergy, <:Any}
}

@inline function Jutul.face_flux!(Q, left, right, face, face_sign, eq::EnergyConservationLaw, state, model, dt, flow_disc::TwoPointPotentialFlowHardCoded)
    # Specific version for tpfa flux
    # TODO: Add general version for thermal
    grad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    # TODO: This could be inconsistent if we use flux_type for something.
    flux_type = Jutul.flux_type(eq)
    q = energy_flux(face, state, model, grad, upw, flux_type, eq)
    return setindex(Q, q, 1)
end

@inline function Jutul.face_flux!(q_i, face, eq::EnergyConservationLaw, state, model, dt, flow_disc::PotentialFlow, ldisc)
    # Inner version, for generic flux
    kgrad, upw = ldisc.face_disc(face)
    ft = Jutul.flux_type(eq)
    q = energy_flux(face, state, model, kgrad, upw, ft, eq)
    return setindex(q_i, q, 1)
end

function energy_flux(face, state, model, grad, upw, flux_type, eq::ConservationLaw{:TotalThermalEnergy})
    return thermal_heat_flux(face, state, model, grad, upw, flux_type)
end

function energy_flux(face, state, model, grad, upw, flux_type, eq::ConservationLaw{:TotalEnergy})
    q = thermal_heat_flux(face, state, model, grad, upw, flux_type)
    return q + potential_energy_flux(face, state, model, grad, upw, flux_type)
end

"""
    potential_energy_flux(face, state, model, grad, upw, flux_type)

Advective flux of gravitational potential energy over a face, i.e. the sum
over phases of the phase mass flux times the upwinded potential energy per unit
mass (see [`UnitPotentialEnergy`](@ref)).
"""
function potential_energy_flux(face, state, model, grad, upw, flux_type)
    Φ = state.UnitPotentialEnergy
    mass_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, grad, upw)
    q = zero(eltype(mass_fluxes))
    for F_α in mass_fluxes
        q += F_α*upwind(upw, Φ, F_α)
    end
    return q
end

"""
    thermal_heat_flux(face, state, model, grad, upw, flux_type)

Calculate the thermal heat flux for a given face in a thermal model.

# Arguments
- `face`: The face for which the heat flux is being calculated.
- `state`: The current state of the system.
- `model`: The thermal model being used.
- `grad`: The gradient operator.
- `upw`: Upwind scheme operator.
- `flux_type`: The type of flux calculation to be used.

# Returns
- The calculated thermal heat flux for the given face.
"""
function thermal_heat_flux(face, state, model, grad, upw, flux_type)
    T = state.Temperature
    H_f = state.FluidEnthalpy
    λ_r = state.RockThermalConductivities
    λ_f = state.FluidThermalConductivities
    S = state.Saturations
    nph = number_of_phases(model.system)

    convective_flux = 0.0
    mass_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, grad, upw)
    for α in 1:nph
        F_α = mass_fluxes[α]
        H_face_α = phase_upwind(upw, H_f, α, F_α)
        convective_flux += H_face_α*F_α
    end

    λ_total = λ_r[face]
    for α in 1:nph
        λ_total += λ_f[α, face]*phase_face_average(S, grad, α)
    end
    conductive_flux = -λ_total*gradient(T, grad)
    return conductive_flux + convective_flux
end


"""
    Jutul.convergence_criterion(model, storage, eq::ConservationLaw{:TotalThermalEnergy}, eq_s, r; dt = 1.0, update_report = missing)

Calculate the convergence criterion for the total thermal energy conservation
law. The same criterion is used for the `:TotalEnergy` formulation.

# Arguments
- `model`: The model object containing the simulation parameters and state.
- `storage`: The storage object used to keep track of intermediate results.
- `eq::ConservationLaw{:TotalThermalEnergy}`: The conservation law for total thermal energy.
- `eq_s`: The storage of the conservation law equation.
- `r`: The residual of the conservation law equation.
- `dt`: The time step size in seconds (default is 1.0).
- `update_report`: An optional argument for updating the report (default is `missing`).

# Returns
- The convergence criterion values for the total thermal energy conservation law (maximum).
"""
function Jutul.convergence_criterion(model, storage, eq::EnergyConservationLaw, eq_s, r; dt = 1.0, update_report = missing)
    M = global_map(model.domain)
    # Note: The total energy formulation is also scaled by the thermal energy,
    # as the potential energy depends on the choice of datum and can be zero or
    # negative.
    E0 = Jutul.active_view(storage.state0.TotalThermalEnergy, M;
        for_variables = false)
    residual = view(r, 1, :)
    Etot = sum(value, E0)
    eb = abs(sum(value, residual))*dt/Etot
    normalized_energy(residual, total_energy) =
        dt*abs(value(residual))/value(total_energy)
    cnv = mapreduce(normalized_energy, max, residual, E0)
    ΔT = temperature_increment(model, storage.state, update_report)

    return (
        CNV = (errors = (cnv, ), names = ("Max", )),
        EB = (errors = (eb, ), names = ("Energy balance", )), 
        increment_dT = (errors = (ΔT, ), names = ("ΔT", )),
        )
end

function temperature_increment(model, state, update_report)
    t_report = update_report[:Temperature]
    if haskey(t_report, :max)
        v = t_report.max
    else
        v = 1.0
    end
    return v
end

function temperature_increment(model, state, update_report::Missing)
    return 1.0
end
