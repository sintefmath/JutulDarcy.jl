
abstract type AbstractReservoirFromWellCT <: Jutul.AdditiveCrossTerm end
struct ReservoirFromWellFlowCT{T<:AbstractVector, I<:AbstractVector} <: AbstractReservoirFromWellCT
    WI::T
    reservoir_cells::I
    well_cells::I
end

function Base.show(io::IO, d::ReservoirFromWellFlowCT)
    n = length(d.WI)
    print(io, "ReservoirFromWellFlowCT ($n connections)")
end

Jutul.symmetry(::AbstractReservoirFromWellCT) = Jutul.CTSkewSymmetry()
Jutul.can_impact_cross_term(force_t::PerforationMask, cross_term::AbstractReservoirFromWellCT) = true

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::ReservoirFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i))
    # Unpack properties
    sys = flow_system(model_t.system)
    @inbounds begin 
        reservoir_cell = ct.reservoir_cells[i]
        well_cell = ct.well_cells[i]
        WI = state_s.WellIndices[i]
        gdz = state_s.PerforationGravityDifference[i]
        p_well = state_s.Pressure
        p_res = state_t.Pressure
        dp = p_well[well_cell] - p_res[reservoir_cell]
    end
    rhoS = reference_densities(sys)

    # Wrap the key connection data in tuple for easy extension later
    conn = (
        dp = dp,
        WI = WI,
        gdz = gdz,
        well = well_cell,
        perforation = i,
        reservoir = reservoir_cell
    )
    # Call smaller interface that is easy to specialize
    if haskey(state_s, :MassFractions)
        @inbounds simple_well_perforation_flux!(out, sys, state_t, state_s, rhoS, conn)
    else
        @inbounds multisegment_well_perforation_flux!(out, sys, state_t, state_s, rhoS, conn)
    end
end

function perforation_phase_potential_difference(conn, state_res, state_well, ix)
    dp = conn.dp
    WI = conn.WI
    if haskey(state_well, :ConnectionPressureDrop)
        dp += state_well.ConnectionPressureDrop[conn.perforation]
    elseif conn.gdz != 0
        ρ_r = state_res.PhaseMassDensities[ix, conn.reservoir]
        if haskey(state_well, :PhaseMassDensities)
            ρ_w = state_well.PhaseMassDensities[ix, conn.well]
            ρ = 0.5*(ρ_r + ρ_w)
        else
            ρ = ρ_r
        end
        dp += ρ*conn.gdz
    end
    return -WI*dp
end

function Jutul.cross_term_entities(ct::AbstractReservoirFromWellCT, eq::ConservationLaw, model)
    return ct.reservoir_cells
end

function Jutul.cross_term_entities_source(ct::AbstractReservoirFromWellCT, eq::ConservationLaw, model)
    return ct.well_cells
end

function Jutul.subcrossterm(ct::ReservoirFromWellFlowCT, ctp, m_t, m_s, map_res::FiniteVolumeGlobalMap, ::TrivialGlobalMap, partition)
    (; WI, reservoir_cells, well_cells) = ct
    rc = map(
        c -> Jutul.local_cell(c, map_res),
        reservoir_cells)
    return ReservoirFromWellFlowCT(copy(WI), rc, copy(well_cells))
end

# Well influence on facility
struct FacilityFromWellFlowCT <: Jutul.AdditiveCrossTerm
    well::Symbol
end

well_top_node() = 1

Jutul.cross_term_entities(ct::FacilityFromWellFlowCT, eq::ControlEquationWell, model) = get_well_position(model.domain, ct.well)

import Jutul: prepare_cross_term_in_entity!

function Jutul.prepare_cross_term_in_entity!(i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i))
    # Check the limits before we calculate the cross term. Then, we know the current control
    # is within limits when it is time to update the cross term itself.
    well_symbol = ct.well
    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    target = ctrl.target
    if !isa(target, DisabledTarget)
        limits = current_limits(cfg, well_symbol)
        if !isnothing(limits)
            rhoS, S = surface_density_and_volume_fractions(state_well)
            q_t = facility_surface_mass_rate_for_well(facility, well_symbol, state_facility)
            apply_well_limit!(cfg, target, well, state_well, well_symbol, rhoS, S, value(q_t), limits)
        end
    end
end

function Jutul.apply_force_to_cross_term!(ct_s, cross_term::ReservoirFromWellFlowCT, target, source, model, storage, dt, force::PerforationMask; time = time)
    mask = force.values
    apply_perforation_mask!(ct_s.target, mask)
    apply_perforation_mask!(ct_s.source, mask)
end

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i))

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
    rhoS, S = surface_density_and_volume_fractions(state_well)
    t = well_target(ctrl, target, well, state_well, rhoS, S)
    if rate_weighted(target)
        actual_rate = t*q_t
        name = physical_representation(well.domain).name
        if abs(actual_rate) < MIN_ACTIVE_WELL_RATE
            t = q_t
        else
            t = actual_rate
        end
    end
    t += 1e-20*q_t
    t_num = target.value
    return (t, t_num)
end

# Facility influence on well
struct WellFromFacilityFlowCT <: Jutul.AdditiveCrossTerm
    well::Symbol
end

Jutul.cross_term_entities(ct::WellFromFacilityFlowCT, eq::ConservationLaw, model) = [well_top_node()]

function update_cross_term_in_entity!(out, i,
    state_well, state0_well,
    state_facility, state0_facility,
    well, facility,
    ct::WellFromFacilityFlowCT, eq, dt, ldisc = local_discretization(ct, i))

    well_symbol = ct.well
    pos = get_well_position(facility.domain, well_symbol)

    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    qT = state_facility.TotalSurfaceMassRate[pos] 
    # Hack for sparsity detection
    qT += 0*bottom_hole_pressure(state_well)
    if ctrl isa DisabledControl
        factor = 1.0
    else
        factor = ctrl.factor
    end
    qT *= factor
    if isa(ctrl, InjectorControl)
        if value(qT) < 0
            @warn "Injector $well_symbol is producing?"
        end
        mix = ctrl.injection_mixture
        nmix = length(mix)
        ncomp = number_of_components(flow_system(well.system))
        @assert nmix == ncomp "Injection composition length ($nmix) must match number of components ($ncomp)."
    else
        if value(qT) > 0 && ctrl isa ProducerControl
            @warn "Producer $well_symbol is injecting?"
        end
        if haskey(state_well, :MassFractions)
            mix = state_well.MassFractions
        else
            masses = @views state_well.TotalMasses[:, well_top_node()]
            mass = sum(masses)
            mix = masses./mass
        end
    end
    for i in eachindex(out)
        @inbounds out[i] = -mix[i]*qT
    end
end

# Thermal
struct ReservoirFromWellThermalCT{T<:AbstractVector, I<:AbstractVector} <: AbstractReservoirFromWellCT
    CI::T
    WI::T
    reservoir_cells::I
    well_cells::I
end

function update_cross_term_in_entity!(out, i,
    state_res, state0_res,
    state_well, state0_well, 
    model_res, model_well,
    ct::ReservoirFromWellThermalCT, eq, dt, ldisc = local_discretization(ct, i))
    # Unpack properties
    sys = thermal_system(model_res.system)
    nph = number_of_phases(sys)
    @inbounds begin 
        reservoir_cell = ct.reservoir_cells[i]
        well_cell = ct.well_cells[i]
        CI = ct.CI[i]
        WI = state_well.WellIndices[i]
        gdz = state_well.PerforationGravityDifference[i]
        p_well = state_well.Pressure
        p_res = state_res.Pressure
        dp = p_well[well_cell] - p_res[reservoir_cell]
        conn = (
            dp = dp,
            WI = WI,
            gdz = gdz,
            well = well_cell,
            perforation = i,
            reservoir = reservoir_cell
        )
    end

    # p_well = state_well.Pressure[well_cell]
    # p_res = state_res.Pressure[reservoir_cell]

    kr = state_res.RelativePermeabilities
    mu = state_res.PhaseViscosities

    λ_t = 0
    for ph in 1:nph
        λ_t += state_res.PhaseMobilities[ph, reservoir_cell]
    end
    advective_heat_flux = 0
    for ph in 1:nph
        q_ph = perforation_phase_mass_flux(λ_t, conn, state_res, state_well, ph)
        if q_ph < 0
            # Injection
            H_perf = state_well.FluidEnthalpy[ph, well_cell]
        else
            H_perf = state_res.FluidEnthalpy[ph, reservoir_cell]
        end
        advective_heat_flux += H_perf*q_ph
    end
    T_well = state_well.Temperature[well_cell]
    T_res = state_res.Temperature[reservoir_cell]

    conductive_heat_flux = -CI*(T_well - T_res)
    out[] = advective_heat_flux + conductive_heat_flux
end

function Base.show(io::IO, d::ReservoirFromWellThermalCT)
    n = length(d.CI)
    print(io, "ReservoirFromWellThermalCT ($n connections)")
end

function Jutul.subcrossterm(ct::ReservoirFromWellThermalCT, ctp, m_t, m_s, map_res::FiniteVolumeGlobalMap, ::TrivialGlobalMap, partition)
    (; CI, WI, reservoir_cells, well_cells) = ct
    rc = map(
        c -> Jutul.local_cell(c, map_res),
        reservoir_cells)
    return ReservoirFromWellThermalCT(copy(CI), copy(WI), rc, copy(well_cells))
end

struct WellFromFacilityThermalCT <: Jutul.AdditiveCrossTerm
    well::Symbol
end

Jutul.cross_term_entities(ct::WellFromFacilityThermalCT, eq::ConservationLaw, model) = [well_top_node()]

function update_cross_term_in_entity!(out, i,
    state_well, state0_well,
    state_facility, state0_facility,
    well, facility,
    ct::WellFromFacilityThermalCT, eq, dt, ldisc = local_discretization(ct, i))
    well_symbol = ct.well
    pos = get_well_position(facility.domain, well_symbol)

    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    qT = state_facility.TotalSurfaceMassRate[pos] 
    # Hack for sparsity detection
    qT += 0*bottom_hole_pressure(state_well)

    cell = well_top_node()
    H = well_top_node_enthalpy(ctrl, state_well, cell)
    out[] = -qT*H
end

function well_top_node_enthalpy(ctrl::InjectorControl, state_well, cell)
    heat_capacity = state_well.FluidHeatCapacity[cell]
    p = state_well.Pressure[cell]
    density = ctrl.mixture_density
    T = ctrl.temperature
    nph = size(state_well.Saturations, 1)
    H = 0
    for ph in 1:nph
        C = state_well.FluidHeatCapacity[ph, cell]
        dens = state_well.PhaseMassDensities[ph, cell]
        S = state_well.Saturations[ph, cell]
        H += S*(C*T + p/dens)
    end
    return H
end

function well_top_node_enthalpy(ctrl, state_well, cell)
    H = state_well.FluidEnthalpy
    S = state_well.Saturations
    H_w = 0.0
    for ph in axes(H, 1)
        H_w += H[ph, cell]*S[ph, cell]
    end
    return H_w
end
