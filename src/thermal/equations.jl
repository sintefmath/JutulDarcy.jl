@inline function Jutul.face_flux!(Q, left, right, face, face_sign, eq::ConservationLaw{:TotalThermalEnergy, <:Any}, state, model, dt, flow_disc::TwoPointPotentialFlowHardCoded)
    # Specific version for tpfa flux
    # TODO: Add general version for thermal
    grad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    # TODO: This could be inconsistent if we use flux_type for something.
    flux_type = Jutul.flux_type(eq)
    q = thermal_heat_flux(face, state, model, grad, upw, flux_type)
    return setindex(Q, q, 1)
end

@inline function Jutul.face_flux!(q_i, face, eq::ConservationLaw{:TotalThermalEnergy, <:Any}, state, model, dt, flow_disc::PotentialFlow, ldisc)
    # Inner version, for generic flux
    kgrad, upw = ldisc.face_disc(face)
    ft = Jutul.flux_type(eq)
    q = thermal_heat_flux(face, state, model, kgrad, upw, ft)
    return setindex(q_i, q, 1)
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
    θ_r = state.RockThermalTransmissibilites
    S = state.Saturations
    nph = number_of_phases(model.system)

    convective_flux = 0.0
    mass_fluxes = darcy_phase_mass_fluxes(face, state, model, flux_type, grad, upw)
    for α in 1:nph
        F_α = mass_fluxes[α]
        H_face_α = phase_upwind(upw, H_f, α, F_α)
        convective_flux += H_face_α*F_α
    end

    # Compute fluid thermal transmissibilities dynamically when conductivity is available.
    has_dynamic_λf = haskey(model.secondary_variables, :FluidThermalConductivity) && haskey(state, :FluidThermalConductivity)
    if has_dynamic_λf
        λ_f = state.FluidThermalConductivity
        θ_geo = state.HalfFaceTransGeo
    else
        θ_f = state.FluidThermalTransmissibilites
    end

    θ = θ_r[face]
    for α in 1:nph
        if has_dynamic_λf
            grad isa TPFA || throw(ArgumentError("TPFA gradient expected for dynamic fluid thermal conductivity"))
            tpfa = grad
            left = tpfa.left
            right = tpfa.right
            ϕ_l = state.FluidVolume[left]./state.BulkVolume[left]
            ϕ_r = state.FluidVolume[right]./state.BulkVolume[right]
            λ_l = ϕ_l.*fluid_phase_cell_value(λ_f, α, left)
            λ_r = ϕ_r.*fluid_phase_cell_value(λ_f, α, right)
            θ_l = @inbounds θ_geo[1, face]*λ_l
            θ_r = @inbounds θ_geo[2, face]*λ_r
            if θ_l == 0.0 || θ_r == 0.0
                θ_fα = 0.0
            else
                θ_fα = 1/(1/θ_l + 1/θ_r)
            end

            # λ_face_α = @inbounds λ_f[α, face]
        else
            θ_fα = θ_f[α, face]
        end
        θ += θ_fα*phase_face_average(S, grad, α)
    end

    if haskey(state, :ThermalDispersionTransmissibility)
        θ_d = state.ThermalDispersionTransmissibility
        @inbounds θ_dl, θ_dr = θ_d[1, face], θ_d[2, face]
        if θ_dl != 0.0 && θ_dr != 0.0
            # Dispersive conductivity is ρ c_p D for the (single) fluid phase
            l, r = grad.left, grad.right
            ρCp_l = @inbounds state.PhaseMassDensities[1, l]*state.ComponentHeatCapacity[1, l]
            ρCp_r = @inbounds state.PhaseMassDensities[1, r]*state.ComponentHeatCapacity[1, r]
            θ += 1/(1/(θ_dl*ρCp_l) + 1/(θ_dr*ρCp_r))
        end
    end

    conductive_flux = -θ*gradient(T, grad)
    return conductive_flux + convective_flux
end

function fluid_phase_cell_value(λ_f, α, cell)
    if λ_f isa AbstractVector
        return @inbounds λ_f[cell]
    elseif λ_f isa AbstractArray
        return @inbounds λ_f[α, cell]
    end
end

"""
    Jutul.convergence_criterion(model, storage, eq::ConservationLaw{:TotalThermalEnergy}, eq_s, r; dt = 1.0, update_report = missing)

Calculate the convergence criterion for the total thermal energy conservation law.

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
function Jutul.convergence_criterion(model, storage, eq::ConservationLaw{:TotalThermalEnergy}, eq_s, r; dt = 1.0, update_report = missing)
    M = global_map(model.domain)
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
