struct PressureReservoirFromWellFlowCT{T} <: Jutul.AdditiveCrossTerm
    parent::T
end

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t,
        state_s, state0_s, 
        model_t, model_s,
        ct::PressureReservoirFromWellFlowCT,
        eq::PressureEquation, dt, ldisc = local_discretization(ct, i)
    )
    # Target is reservoir, source is well
    sys = model_t.system
    rhoS = reference_densities(sys)
    conn = cross_term_perforation_get_conn(ct.parent, i, state_s, state_t)

    ncomp = number_of_components(sys)
    out_full = zeros(eltype(out), ncomp)
    @assert !haskey(state_s, :MassFractions)
    q = multisegment_well_perforation_flux!(out_full, sys, state_t, state_s, rhoS, conn)
    val = zero(eltype(out))
    for i in eachindex(q)
        val += q[i]*state_t.PressureReductionFactors[i, conn.reservoir]
    end
    out[1] = val
end

struct PressureWellFromReservoirFlowCT{T} <: Jutul.AdditiveCrossTerm
    parent::T
end

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t,
        state_s, state0_s, 
        model_t, model_s,
        ct::PressureWellFromReservoirFlowCT,
        eq::PressureEquation, dt, ldisc = local_discretization(ct, i)
    )
    # Target is well, source is reservoir
    sys = model_t.system
    rhoS = reference_densities(sys)
    conn = cross_term_perforation_get_conn(ct.parent, i, state_t, state_s)

    ncomp = number_of_components(sys)
    out_full = zeros(eltype(out), ncomp)
    q = multisegment_well_perforation_flux!(out_full, sys, state_s, state_t, rhoS, conn)
    val = zero(eltype(out))
    for i in eachindex(q)
        val += q[i]*state_t.PressureReductionFactors[i, conn.well]
    end
    out[1] = -val
end

function Jutul.cross_term_entities(ct::PressureReservoirFromWellFlowCT, eq::PressureEquation, model)
    return ct.parent.reservoir_cells
end

function Jutul.cross_term_entities_source(ct::PressureReservoirFromWellFlowCT, eq::PressureEquation, model)
    return ct.parent.well_cells
end

function Jutul.cross_term_entities(ct::PressureWellFromReservoirFlowCT, eq::PressureEquation, model)
    return ct.parent.well_cells
end

function Jutul.cross_term_entities_source(ct::PressureWellFromReservoirFlowCT, eq::PressureEquation, model)
    return ct.parent.reservoir_cells
end
