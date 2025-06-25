struct PressureReservoirFromWellFlowCT{T} <: Jutul.AdditiveCrossTerm
    parent::T
end

function Jutul.update_cross_term_in_entity!(out, i,
        state_res, state0_res,
        state_well, state0_well, 
        model_res, model_well,
        ct::PressureReservoirFromWellFlowCT,
        eq::PressureEquation, dt, ldisc = local_discretization(ct, i)
    )
    # Target is reservoir, source is well
    model_well::SimpleWellModel
    sys = model_res.system
    dest = state_res
    T = eltype(out)
    out[1] = pressure_perforation_flux(T, ct.parent, i, dest, state_res, state_well, sys, false)
end

struct PressureWellFromReservoirFlowCT{T} <: Jutul.AdditiveCrossTerm
    parent::T
end

function Jutul.update_cross_term_in_entity!(out, i,
        state_well, state0_t,
        state_res, state0_s, 
        model_well, model_res,
        ct::PressureWellFromReservoirFlowCT,
        eq::PressureEquation, dt, ldisc = local_discretization(ct, i)
    )
    # Target is well, source is reservoir
    model_well::SimpleWellModel
    sys = model_well.system
    dest = state_well
    T = eltype(out)
    out[1] = -pressure_perforation_flux(T, ct.parent, i, dest, state_res, state_well, sys, true)
end

function pressure_perforation_flux(T, ct, i, state_dest, state_res, state_well, sys, is_well)
    conn = cross_term_perforation_get_conn(ct, i, state_well, state_res)
    rhoS = reference_densities(sys)
    ncomp = number_of_components(sys)
    q = zeros(T, ncomp)
    multisegment_well_perforation_flux!(q, sys, state_res, state_well, rhoS, conn)
    val = zero(T)
    w = state_dest.PressureReductionFactors
    @assert length(q) == size(w, 1) "Size mismatch: expected $(size(w, 1)), got $(length(q)) for $ncomp components."
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
