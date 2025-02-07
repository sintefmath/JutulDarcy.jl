function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel{<:Any, T}, eq::ConservationLaw{:TotalMasses}, eq_s, force::V, time) where {T<:MultiPhaseCompositionalSystemLV, V <: AbstractVector{<:FlowBoundaryCondition}}
    state = storage.state
    p = state.Pressure
    for bc in force
        c = bc.cell
        acc_i = view(acc, :, c)
        T_f = bc.trans_flow
        Δp = p[c] - bc.pressure
        q = T_f*Δp
        apply_flow_bc!(acc_i, q, bc, model, state, time)
    end
end


function apply_flow_bc!(acc, q, bc, model::SimulationModel{<:Any, T}, state, time) where T<:MultiPhaseCompositionalSystemLV
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
        q_v = rho[l, c]*mob[l, c]*q
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
                acc[c] += q*rho_inj*λ_t*F
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

