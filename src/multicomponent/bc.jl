# function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel{<:Any, T}, eq::ConservationLaw{:TotalMasses}, eq_s, force::V, time) where {T<:AbstractCompositionalSystemLV, V <: AbstractVector{<:FlowBoundaryCondition}}
#     state = storage.state
#     p = state.Pressure
#     for bc in force
#         c = bc.cell
#         acc_i = view(acc, :, c)
#         q = compute_bc_mass_fluxes(model.system, bc, global_map(model), state)
#         apply_flow_bc!(acc_i, q, bc, model, state, time)
#     end
# end

function compute_bc_mass_fluxes(system::AbstractCompositionalSystemLV, bc, gmap, state)
    c = bc.cell
    T_f = bc.trans_flow
    p = state.Pressure
    Δp = p[c] - bc.pressure
    q_vol = T_f*Δp

    mob = state.PhaseMobilities
    rho = state.PhaseMassDensities
    ncomp = number_of_components(system)
    nph = number_of_phases(system)

    rho_inj = bc.density
    f_inj = bc.fractional_flow
    c = bc.cell
    # TODO: Capillary pressure.
    T = Base.promote_type(typeof(q_vol), eltype(rho), eltype(mob))
    if q_vol > 0
        q = zeros(T, nph, ncomp)
        X = state.LiquidMassFractions
        Y = state.VaporMassFractions
        # Pressure inside is higher than outside, flow out from domain
        phase_ix = phase_indices(system)
        if has_other_phase(system)
            a, l, v = phase_ix
            ncomp_mix = ncomp-1
            q[a, ncomp_mix+1] = rho[a, c]*mob[a, c]*q_vol
        else
            ncomp_mix = ncomp
            l, v = phase_ix
        end
        q_l = rho[l, c]*mob[l, c]*q_vol
        q_v = rho[v, c]*mob[v, c]*q_vol
        for i in 1:ncomp_mix
            q[l, i] = q_l*X[i, c]
            q[v, i] = q_v*Y[i, c]
        end
    else
        q = zeros(T, 1, ncomp)
        # TODO: This is duplicated code, factor out...
        # Injection of mass
        λ_t = 0.0
        for ph in 1:nph
            λ_t += mob[ph, c]
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
            for i in 1:ncomp
                total += state.TotalMasses[i, c]
            end
            for i in 1:ncomp
                F = state.TotalMasses[i, c]/total
                q[i] = q_vol*rho_inj*λ_t*F
            end
        else
            @assert length(f_inj) == ncomp
            for i in 1:ncomp
                F = f_inj[i]
                q[i] = q_vol*rho_inj*λ_t*F
            end
        end
    end
    return q

end

function compute_bc_heat_fluxes(system::AbstractCompositionalSystemLV, bc, gmap, state, nph)

    nph = number_of_phases(system)
    q = compute_bc_mass_fluxes(system, bc, gmap, state)
    c = Jutul.full_cell(bc.cell, gmap)

    # Get reservoir properties
    h_ph = state.FluidEnthalpy
    h = state.Enthalpy

    q_adv = 0.0
    if size(q, 1) == 1
        q_tot = sum(q)
        @assert q_tot <= 0.0
        q_adv = q_tot*bc.enthalpy
    else
        @assert size(q, 1) == nph
        for ph in 1:nph
            q_ph = sum(q[ph, :])
            q_adv += q_ph*h_ph[ph, c]
        end
    end

    T_h    = bc.trans_thermal
    T      = state.Temperature
    T_bc   = bc.temperature
    ΔT = T[c] - T_bc
    qh_cond = T_h*ΔT

    return q_adv, qh_cond

end

function apply_flow_bc!(acc, q, bc, model::SimulationModel{<:Any, T}, state, time) where T<:AbstractCompositionalSystemLV

    system = reservoir_model(model).system
    ncomp = number_of_components(system)
    for i in 1:ncomp
        acc[i] += sum(q[:, i])
    end

end

function apply_flow_bc!(acc, q, bc, model::SimulationModel{<:Any, T}, state, time) where T<:AbstractCompositionalSystemLV
    mob = state.PhaseMobilities
    rho = state.PhaseMassDensities
    ncomp = length(acc)
    nph = size(rho, 1)

    rho_inj = bc.density
    f_inj = bc.fractional_flow
    c = bc.cell
    # TODO: Capillary pressure.
    if q > 0
        X = state.LiquidMassFractions
        Y = state.VaporMassFractions
        # Pressure inside is higher than outside, flow out from domain
        sys = model.system
        phase_ix = phase_indices(sys)
        if has_other_phase(sys)
            a, l, v = phase_ix
            ncomp_mix = ncomp-1
            acc[ncomp_mix+1] += rho[a, c]*mob[a, c]*q
        else
            ncomp_mix = ncomp
            l, v = phase_ix
        end
        q_l = rho[l, c]*mob[l, c]*q
        q_v = rho[v, c]*mob[v, c]*q
        for i in 1:ncomp_mix
            acc[i] += q_l*X[i, c] + q_v*Y[i, c]
        end
    else
        # TODO: This is duplicated code, factor out...
        # Injection of mass
        λ_t = 0.0
        for ph in 1:nph
            λ_t += mob[ph, c]
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
            for i in 1:ncomp
                total += state.TotalMasses[i, c]
            end
            for i in 1:ncomp
                F = state.TotalMasses[i, c]/total
                acc[i] += q*rho_inj*λ_t*F
            end
        else
            @assert length(f_inj) == ncomp
            for i in 1:ncomp
                F = f_inj[i]
                acc[i] += q*rho_inj*λ_t*F
            end
        end
    end
end

function compute_bc_heat_fluxes(bc, system, gmap, state)

    # function apply_flow_bc!(acc, q, bc, model::SimulationModel{<:Any, T}, state, time) where T<:AbstractCompositionalSystemLV
    mob = state.PhaseMobilities
    rho = state.PhaseMassDensities
    h_ph = state.FluidEnthalpy
    h = state.Enthalpy
    ncomp = number_of_components(system)
    nph = size(rho, 1)

    rho_inj = bc.density
    f_inj = bc.fractional_flow
    c = bc.cell
    # TODO: Capillary pressure.
    if q > 0
        X = state.LiquidMassFractions
        Y = state.VaporMassFractions
        # Pressure inside is higher than outside, flow out from domain
        sys = model.system
        phase_ix = phase_indices(sys)
        if has_other_phase(sys)
            a, l, v = phase_ix
            ncomp_mix = ncomp-1
            acc[ncomp_mix+1] += rho[a, c]*mob[a, c]*q
        else
            ncomp_mix = ncomp
            l, v = phase_ix
        end
        q_l = rho[l, c]*mob[l, c]*q
        q_v = rho[v, c]*mob[v, c]*q
        q_adv = q_l*h_ph[l, c] + q_v*h_ph[v, c]
        # for i in 1:ncomp_mix
        #     acc[i] += h_ph[l, c]*q_l*X[i, c] + h_ph[v, c]*q_v*Y[i, c]
        # end
    else
        # TODO: This is duplicated code, factor out...
        # Injection of mass
        λ_t = 0.0
        for ph in 1:nph
            λ_t += mob[ph, c]
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
            for i in 1:ncomp
                total += state.TotalMasses[i, c]
            end

            for i in 1:ncomp
                F = state.TotalMasses[i, c]/total
                acc[i] += h[c]q*rho_inj*λ_t*F
            end
        else
            @assert length(f_inj) == ncomp
            for i in 1:ncomp
                F = f_inj[i]
                acc[i] += h[c]*q*rho_inj*λ_t*F
            end
        end
    end
end