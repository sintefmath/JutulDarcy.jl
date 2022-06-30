struct ReservoirWellCrossTerm{T<:AbstractVector, I<:AbstractVector} <: Jutul.AdditiveCrossTerm
    WI::T
    reservoir_cells::I
    well_cells::I
end

Jutul.symmetry(::ReservoirWellCrossTerm) = Jutul.CTSkewSymmetry()

function update_cross_term_in_entity!(out, i, state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    param_t, param_s,
    ct::ReservoirWellCrossTerm, eq, dt, ldisc = local_discretization(ct, i))
    # Unpack properties
    reservoir_cell = ct.reservoir_cells[i]
    well_cell = ct.well_cells[i]
    WI = ct.WI[i]
    rhoS = param_s[:reference_densities]
    # Call smaller interface that is easy to specialize
    well_perforation_flux!(out, model_t.system, state_t, state_s, rhoS, WI, reservoir_cell, well_cell)
end

Jutul.cross_term_entities(ct::ReservoirWellCrossTerm, eq::ConservationLaw, model) = ct.reservoir_cells
Jutul.cross_term_entities_source(ct::ReservoirWellCrossTerm, eq::ConservationLaw, model) = ct.well_cells
