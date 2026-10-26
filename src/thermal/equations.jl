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
    θ_r = state.RockThermalConductivities
    S = state.Saturations
    nph = number_of_phases(model.system)

    convective_flux = 0.0
    mass_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, grad, upw)
    for α in 1:nph
        F_α = mass_fluxes[α]
        H_face_α = phase_upwind(upw, H_f, α, F_α)
        convective_flux += H_face_α*F_α
    end

    θ_f = state.FluidThermalConductivities
    θ = θ_r[face]
    for α in 1:nph
        θ_fα = fluid_thermal_transmissibility(state, θ_f, face, grad, α)
        θ += θ_fα*phase_face_average(S, grad, α)
    end
    conductive_flux = -θ*gradient(T, grad)
    return conductive_flux + convective_flux
end

"""
    tpfa_face_cells(tpfa::TPFA)

Cells of a two-point gradient ordered as in the face neighborship, i.e.
`(N[1, face], N[2, face])`, independent of which cell the flux is evaluated
for.
"""
@inline function tpfa_face_cells(tpfa::TPFA)
    if tpfa.face_sign > 0
        return (tpfa.left, tpfa.right)
    else
        return (tpfa.right, tpfa.left)
    end
end

"""
    fluid_thermal_transmissibility(state, θ_f, face, grad, α)

Fluid thermal transmissibility of phase `α` over `face`. If the model has a
cell-wise fluid thermal conductivity (see
`set_dynamic_fluid_thermal_conductivity!`), the static transmissibility is
scaled by the current conductivity using `FluidThermalHalfFaceWeights`.
"""
@inline function fluid_thermal_transmissibility(state, θ_f, face, grad, α)
    θ_fα = @inbounds θ_f[α, face]
    if haskey(state, :FluidThermalHalfFaceWeights)
        θ_fα *= fluid_thermal_conductivity_factor(state.FluidThermalConductivity, state.FluidThermalHalfFaceWeights, face, grad, α)
    end
    return θ_fα
end

@inline function fluid_thermal_conductivity_factor(λ, w, face, tpfa::TPFA, α)
    l, r = tpfa_face_cells(tpfa)
    a_l = @inbounds w[2α - 1, face]*λ[α, l]
    a_r = @inbounds w[2α, face]*λ[α, r]
    if a_l > 0 && a_r > 0
        return a_l*a_r/(a_l + a_r)
    else
        return zero(a_l*a_r)
    end
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
