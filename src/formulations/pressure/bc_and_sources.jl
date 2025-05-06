function pressure_accumulation_buffer(acc, eq_s::PressureEquationTPFAStorage)
    return eq_s.buf
end

function pressure_accumulation_buffer(acc, eq_s)
    return acc[:, 1]
end

function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel{D, S}, eq::PressureEquation, eq_s, force::V, time) where {V <: AbstractVector{<:FlowBoundaryCondition}, D, S<:MultiPhaseSystem}
    # Jutul.apply_forces_to_equation!(acc, storage, model, eq.conservation, eq_s, force, time)
    acc = Jutul.get_diagonal_entries(eq, eq_s)
    acc_i = pressure_accumulation_buffer(acc, eq_s)
    state = storage.state
    w = storage.state.PressureReductionFactors
    nph = number_of_phases(reservoir_model(model).system)
    for bc in force
        c = bc.cell
        @. acc_i = 0.0
        q = compute_bc_mass_fluxes(bc, state, nph)
        apply_flow_bc!(acc_i, q, bc, model, state, time)
        for i in eachindex(acc_i)
            acc[i] += w[i, c]*acc_i[i]
        end
    end
    return acc
end

function Jutul.apply_forces_to_equation!(acc, storage, model::SimulationModel{D, S}, eq::PressureEquation, eq_s, force::V, time) where {V <: AbstractVector{SourceTerm{I, F, T}}, D, S<:MultiPhaseSystem} where {I, F, T}
    state = storage.state
    acc = Jutul.get_diagonal_entries(eq, eq_s)
    kr = state.RelativePermeabilities
    mu = state.PhaseViscosities
    w = storage.state.PressureReductionFactors

    rhoS = reference_densities(model.system)
    nph = number_of_phases(model.system)
    M = global_map(model.domain)
    for src in force
        c = Jutul.full_cell(src.cell, M)
        for ph = 1:nph
            q_ph = phase_source(c, src, rhoS[ph], kr, mu, ph)
            @inbounds acc[src.cell] -= w[ph, src.cell]*q_ph
        end
    end
end
