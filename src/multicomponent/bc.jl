function compute_bc_mass_fluxes(system::CompositionalSystemLV, bc, gmap, state)
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

function compute_bc_heat_fluxes(system::CompositionalSystemLV, bc, gmap, state)

    nph = number_of_phases(system)
    q = compute_bc_mass_fluxes(system, bc, gmap, state)
    c = Jutul.full_cell(bc.cell, gmap)

    # Get reservoir properties
    h_ph = state.FluidEnthalpy

    q_adv = 0.0
    if size(q, 1) == 1
        q_tot = sum(q)
        @assert q_tot <= 0.0
        q_adv = q_tot*bc_inflow_enthalpy(bc, state, c)
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

function apply_flow_bc!(acc, q, bc, model::SimulationModel{<:Any, T}, state, time) where T<:CompositionalSystemLV

    system = reservoir_model(model).system
    ncomp = number_of_components(system)
    if ncomp == 1
        acc[] += sum(q)
    else
        for i in 1:ncomp
            acc[i] += sum(q[:, i])
        end
    end

end
