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
    model_s::SimpleWellModel
    sys = model_t.system
    T = eltype(out)
    out[1] = pressure_perforation_flux(T, ct.parent, i, state_t, state_t, state_s, sys, false)
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
    model_t::SimpleWellModel
    sys = model_t.system
    out[1] = -pressure_perforation_flux(eltype(out), ct.parent, i, state_t, state_s, state_t, sys, true)
end

function pressure_perforation_flux(T, ct, i, state_dest, state_res, state_well, sys, is_well)
    conn = cross_term_perforation_get_conn(ct, i, state_well, state_res)
    rhoS = reference_densities(sys)
    ncomp = number_of_components(sys)
    out_full = zeros(T, ncomp)
    q = multisegment_well_perforation_flux!(out_full, sys, state_res, state_well, rhoS, conn)
    val = zero(T)
    w = state_dest.PressureReductionFactors
    @assert length(q) == size(w, 1)
    if is_well
        ix = conn.well
    else
        ix = conn.reservoir
    end
    for i in eachindex(q)
        val += q[i]*w[i, ix]
    end
    return val
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
