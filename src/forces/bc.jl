"""
    FlowBoundaryCondition(
    cell,
    pressure = DEFAULT_MINIMUM_PRESSURE, 
    temperature = 298.15;
    fractional_flow = nothing,
    density = nothing,
    trans_flow = 1e-12,
    trans_thermal = 1e-6
    )

Boundary condition for constant values (pressure/temperature)
"""
function FlowBoundaryCondition(
    cell,
    pressure = DEFAULT_MINIMUM_PRESSURE, 
    temperature = 298.15;
    fractional_flow = nothing,
    density = nothing,
    trans_flow = 1e-12,
    trans_thermal = 1e-6
    )
    @assert isnothing(density) || density > 0.0 "Density, if provided, must be positive"
    if isnothing(fractional_flow)
        f = fractional_flow
    else
        f = Tuple(fractional_flow)
        @assert all(f .> 0)
        @assert sum(f) == 1.0 "Fractional flow for boundary condition in cell $cell must sum to 1."
    end
    @assert pressure >= DEFAULT_MINIMUM_PRESSURE
    @assert temperature >= 0.0
    return FlowBoundaryCondition(cell, pressure, temperature, trans_flow, trans_thermal, f, density)
end

function Jutul.subforce(s::AbstractVector{S}, model) where S<:FlowBoundaryCondition
    error("Not implemented.")
end


function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel{D, S}, eq::ConservationLaw{:TotalMasses}, eq_s, force::V, time) where {V <: AbstractVector{<:FlowBoundaryCondition}, D, S<:MultiPhaseSystem}
    state = storage.state
    p = state.Pressure
    for bc in force
        c = bc.cell
        T_f = bc.trans_flow
        Δp = p[c] - bc.pressure
        q = T_f*Δp
        acc_i = view(acc, :, c)
        apply_flow_bc!(acc_i, q, bc, model, state, time)
    end
end

function apply_flow_bc!(acc, q, bc, model::SimulationModel{<:Any, T}, state, time) where T<:Union{ImmiscibleSystem, SinglePhaseSystem}
    mu = state.PhaseViscosities
    kr = state.RelativePermeabilities
    rho = state.PhaseMassDensities
    nph = length(acc)
    @assert size(kr, 1) == nph

    rho_inj = bc.density
    f_inj = bc.fractional_flow
    c = bc.cell
    if q > 0
        # Pressure inside is higher than outside, flow out from domain
        for ph in eachindex(acc)
            # Immiscible: Density * total flow rate * mobility for each phase
            q_i = q*rho[ph, c]*kr[ph, c]/mu[ph, c]
            acc[ph] += q_i
        end
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
