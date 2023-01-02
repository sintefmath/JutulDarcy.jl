export ReservoirFromWellCT, FacilityFromWellCT, WellFromFacilityCT
struct ReservoirFromWellCT{T<:AbstractVector, I<:AbstractVector} <: Jutul.AdditiveCrossTerm
    WI::T
    reservoir_cells::I
    well_cells::I
end

function Base.show(io::IO, d::ReservoirFromWellCT)
    n = length(d.WI)
    print(io, "ReservoirFromWellCT ($n connections)")
end

Jutul.symmetry(::ReservoirFromWellCT) = Jutul.CTSkewSymmetry()
Jutul.can_impact_cross_term(force_t::PerforationMask, cross_term::ReservoirFromWellCT) = true

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::ReservoirFromWellCT, eq, dt, ldisc = local_discretization(ct, i))
    # Unpack properties
    sys = flow_system(model_t.system)
    @inbounds begin 
        reservoir_cell = ct.reservoir_cells[i]
        well_cell = ct.well_cells[i]
        WI = state_s.WellIndices[i]
    end
    rhoS = reference_densities(sys)

    p_well = state_s.Pressure
    p_res = state_t.Pressure
    # Todo: Fix conn -> cell pressure drop
    ρgdz = 0
    @inbounds dp = -WI*(p_well[well_cell] - p_res[reservoir_cell] + ρgdz)
    # Call smaller interface that is easy to specialize
    @inbounds well_perforation_flux!(out, sys, state_t, state_s, rhoS, dp, reservoir_cell, well_cell)
end

Jutul.cross_term_entities(ct::ReservoirFromWellCT, eq::ConservationLaw, model) = ct.reservoir_cells
Jutul.cross_term_entities_source(ct::ReservoirFromWellCT, eq::ConservationLaw, model) = ct.well_cells

function Jutul.subcrossterm(ct::ReservoirFromWellCT, ctp, m_t, m_s, map_res::FiniteVolumeGlobalMap, ::TrivialGlobalMap, partition)
    (; WI, reservoir_cells, well_cells) = ct
    # rc = map(
    #     c -> Jutul.interior_cell(
    #         Jutul.local_cell(c, map_res),
    #     map_res),
    #     reservoir_cells)
    rc = map(
        c -> Jutul.local_cell(c, map_res),
        reservoir_cells)
    return ReservoirFromWellCT(copy(WI), rc, copy(well_cells))
end

# Well influence on facility
struct FacilityFromWellCT <: Jutul.AdditiveCrossTerm
    well::Symbol
end

well_top_node() = 1

Jutul.cross_term_entities(ct::FacilityFromWellCT, eq::ControlEquationWell, model) = get_well_position(model.domain, ct.well)

function Jutul.prepare_cross_term_in_entity!(i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellCT, eq, dt, ldisc = local_discretization(ct, i))
    # Check the limits before we calculate the cross term. Then, we know the current control
    # is within limits when it is time to update the cross term itself.
    well_symbol = ct.well
    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    target = ctrl.target
    if !isa(target, DisabledTarget)
        limits = current_limits(cfg, well_symbol)
        if !isnothing(limits)
            rhoS = reference_densities(well.system)
            rhoS, S = flash_wellstream_at_surface(well, state_well, rhoS)
            rhoS = tuple(rhoS...)
            q_t = facility_surface_mass_rate_for_well(facility, well_symbol, state_facility)
            apply_well_limit!(cfg, target, well, state_well, well_symbol, rhoS, S, value(q_t), limits)
        end
    end
end

function Jutul.apply_force_to_cross_term!(ct_s, cross_term::ReservoirFromWellCT, target, source, model, storage, dt, force::PerforationMask; time = time)
    mask = force.values
    apply_perforation_mask!(ct_s.target, mask)
    apply_perforation_mask!(ct_s.source, mask)
end

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellCT, eq, dt, ldisc = local_discretization(ct, i))

    well_symbol = ct.well
    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)

    target = ctrl.target
    q_t = facility_surface_mass_rate_for_well(facility, well_symbol, state_facility)
    t, t_num = target_actual_pair(target, well, state_well, q_t, ctrl)
    t += 0*bottom_hole_pressure(state_well) + 0*q_t
    scale = target_scaling(target)
    eq = (t - t_num)/scale
    out[] = eq
end

function target_actual_pair(target::DisabledTarget, well, state_well, q_t, ctrl)
    # The well should have zero rate. Enforce this by the trivial residual R = q_t = 0
    t = q_t
    t_num = 0.0
    return (t, t_num)
end


function target_actual_pair(target, well, state_well, q_t, ctrl)
    need_rates = isa(ctrl, ProducerControl) && !isa(target, BottomHolePressureTarget)
    rhoS = reference_densities(well.system)
    if need_rates
        rhoS, S = flash_wellstream_at_surface(well, state_well, rhoS)
    else
        S = nothing
    end
    rhoS = tuple(rhoS...)
    t = well_target(ctrl, target, well, state_well, rhoS, S)
    if rate_weighted(target)
        t *= q_t
    end
    t_num = target.value
    return (t, t_num)
end

# Facility influence on well
struct WellFromFacilityCT <: Jutul.AdditiveCrossTerm
    well::Symbol
end

Jutul.cross_term_entities(ct::WellFromFacilityCT, eq::ConservationLaw, model) = [well_top_node()]

function update_cross_term_in_entity!(out, i,
    state_well, state0_well,
    state_facility, state0_facility,
    well, facility,
    ct::WellFromFacilityCT, eq, dt, ldisc = local_discretization(ct, i))

    well_symbol = ct.well
    pos = get_well_position(facility.domain, well_symbol)

    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    qT = state_facility.TotalSurfaceMassRate[pos] 
    # Hack for sparsity detection
    qT += 0*bottom_hole_pressure(state_well)

    if isa(ctrl, InjectorControl)
        if value(qT) < 0
            @warn "Injector $well_symbol is producing?"
        end
        mix = ctrl.injection_mixture
        nmix = length(mix)
        ncomp = number_of_components(well.system)
        @assert nmix == ncomp "Injection composition length ($nmix) must match number of components ($ncomp)."
    else
        if value(qT) > 0
            @warn "Producer $well_symbol is injecting?"
        end
        masses = @views state_well.TotalMasses[:, well_top_node()]
        mass = sum(masses)
        mix = masses./mass
    end
    for i in eachindex(out)
        @inbounds out[i] = -mix[i]*qT
    end
end

