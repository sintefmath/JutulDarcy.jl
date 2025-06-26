
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

"""
    update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::ReservoirFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i))

Update mass flow between reservoir and well.
"""
function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    ct::ReservoirFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i))
    sys = model_t.system
    rhoS = reference_densities(sys)
    conn = cross_term_perforation_get_conn(ct, i, state_s, state_t)
    # Call smaller interface that is easy to specialize
    if haskey(state_s, :MassFractions)
        @inbounds simple_well_perforation_flux!(out, sys, state_t, state_s, rhoS, conn)
    else
        @inbounds multisegment_well_perforation_flux!(out, sys, state_t, state_s, rhoS, conn)
    end
end

function cross_term_perforation_get_conn(ct, i, state_s, state_t)
    @inbounds begin 
        reservoir_cell = ct.reservoir_cells[i]
        well_cell = ct.well_cells[i]
        WI = state_s.WellIndices[i]
        gdz = state_s.PerforationGravityDifference[i]
        p_well = state_s.Pressure
        p_res = state_t.Pressure
        dp = p_well[well_cell] - p_res[reservoir_cell]
    end

    # Wrap the key connection data in tuple for easy extension later
    conn = (
        dp = dp,
        WI = WI,
        gdz = gdz,
        well = well_cell,
        perforation = i,
        reservoir = reservoir_cell
    )
    return conn
end

function perforation_phase_potential_difference(conn, state_res, state_well, ix)
    dp = conn.dp
    WI = conn.WI
    WI, dp = Base.promote(WI, dp)
    if haskey(state_res, :PermeabilityMultiplier)
        K_mul = state_res[:PermeabilityMultiplier][conn.reservoir]
        WI *= K_mul
    end
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
            q_t = facility_surface_mass_rate_for_well(facility, well_symbol, state_facility, effective = false)
            apply_well_limit!(cfg, target, well, state_well, well_symbol, rhoS, S, value(q_t), limits)
        end
    end
end

function Jutul.apply_force_to_cross_term!(ct_s, cross_term::ReservoirFromWellFlowCT, target, source, model, storage, dt, force::PerforationMask; time = time)
    mask = force.values
    apply_perforation_mask!(ct_s.target, mask)
    apply_perforation_mask!(ct_s.source, mask)
end

"""
    update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i))

Update the control equation of the facility based on the current well state.
"""
function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellFlowCT, eq, dt, ldisc = local_discretization(ct, i))

    well_symbol = ct.well
    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)

    target = ctrl.target
    update_target!(ctrl, target, state_facility, state_well, facility)
    q_t = facility_surface_mass_rate_for_well(
        facility,
        well_symbol,
        state_facility,
        effective = false
    )
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

"""
    update_cross_term_in_entity!(out, i,
    state_well, state0_well,
    state_facility, state0_facility,
    well, facility,
    ct::WellFromFacilityFlowCT, eq, dt, ldisc = local_discretization(ct, i))

Update the cross-term of the well based on the current facility state. This is
done by adding a source term to the well equation based on the current facility
status (injecting or producing).
"""
function update_cross_term_in_entity!(out, i,
    state_well, state0_well,
    state_facility, state0_facility,
    well, facility,
    ct::WellFromFacilityFlowCT, eq, dt, ldisc = local_discretization(ct, i))

    well_symbol = ct.well
    pos = get_well_position(facility.domain, well_symbol)

    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    q_t = facility_surface_mass_rate_for_well(
        facility,
        well_symbol,
        state_facility,
        effective = true
    )
    # Hack for sparsity detection
    q_t += 0*bottom_hole_pressure(state_well)

    if isa(ctrl, InjectorControl)
        if value(q_t) < 0
            @warn "Injector $well_symbol is producing?"
        end
        mix = ctrl.injection_mixture
        nmix = length(mix)
        ncomp = number_of_components(flow_system(well.system))
        @assert nmix == ncomp "Injection composition length ($nmix) must match number of components ($ncomp)."
    else
        if value(q_t) > 0 && ctrl isa ProducerControl
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
        @inbounds out[i] = -mix[i]*q_t
    end
end

# Thermal
struct ReservoirFromWellThermalCT{T<:AbstractVector, I<:AbstractVector} <: AbstractReservoirFromWellCT
    WIth::T
    WI::T
    reservoir_cells::I
    well_cells::I
end

"""
    update_cross_term_in_entity!(out, i,
    state_res, state0_res,
    state_well, state0_well, 
    model_res, model_well,
    ct::ReservoirFromWellThermalCT, eq, dt, ldisc = local_discretization(ct, i))

Update the cross term between a well and reservoir for thermal equations. This
computes the energy transfer into or out from the well bore and the reservoir,
including both the effect of advection and conduction.
"""
function update_cross_term_in_entity!(out, i,
    state_res, state0_res,
    state_well, state0_well, 
    model_res, model_well,
    ct::ReservoirFromWellThermalCT, eq, dt, ldisc = local_discretization(ct, i))
    # Unpack properties
    sys = model_res.system
    nph = number_of_phases(sys)
    @inbounds begin 
        reservoir_cell = ct.reservoir_cells[i]
        well_cell = ct.well_cells[i]
        WIth = state_well.WellIndicesThermal[i]
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

    kr = state_res.RelativePermeabilities
    mu = state_res.PhaseViscosities

    λ_t = sum(perforation_reservoir_mobilities(state_res, state_well, sys, reservoir_cell, well_cell))
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

    conductive_heat_flux = -WIth*(T_well - T_res)
    out[] = advective_heat_flux + conductive_heat_flux
end

function Base.show(io::IO, d::ReservoirFromWellThermalCT)
    n = length(d.WIth)
    print(io, "ReservoirFromWellThermalCT ($n connections)")
end

function Jutul.subcrossterm(ct::ReservoirFromWellThermalCT, ctp, m_t, m_s, map_res::FiniteVolumeGlobalMap, ::TrivialGlobalMap, partition)
    (; WIth, WI, reservoir_cells, well_cells) = ct
    rc = map(
        c -> Jutul.local_cell(c, map_res),
        reservoir_cells)
    return ReservoirFromWellThermalCT(copy(WIth), copy(WI), rc, copy(well_cells))
end

function Jutul.apply_force_to_cross_term!(ct_s, cross_term::ReservoirFromWellThermalCT, target, source, model, storage, dt, force::PerforationMask; time = time)
    mask = force.values
    apply_perforation_mask!(ct_s.target, mask)
    apply_perforation_mask!(ct_s.source, mask)
end

struct WellFromFacilityThermalCT <: Jutul.AdditiveCrossTerm
    well::Symbol
end

Jutul.cross_term_entities(ct::WellFromFacilityThermalCT, eq::ConservationLaw, model) = [well_top_node()]

"""
    update_cross_term_in_entity!(out, i,
    state_well, state0_well,
    state_facility, state0_facility,
    well, facility,
    ct::WellFromFacilityThermalCT, eq, dt, ldisc = local_discretization(ct, i))

Update the cross term between a well and facility for thermal equations.
"""
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

    T = get_target_temperature(ctrl, ctrl.target, facility, state_facility)
    H = well_top_node_enthalpy(ctrl, well, state_well, T, cell)
    out[] = -qT*H
end

function get_target_temperature(ctrl, target, facility, state_facility)
    return missing
end

function get_target_temperature(ctrl::InjectorControl, target, facility, state_facility)
    return ctrl.temperature
end

function get_target_temperature(ctrl::InjectorControl, target::ReinjectionTarget, facility, state_facility)

    # TODO: This currently assumes constant fluid heat capacity and equal
    # pressures. Should ideally be replaced by enthalpy, which requires
    # FluidEnthalpy to be a Facility variable

    if !isnan(ctrl.temperature)
        return ctrl.temperature
    end

    q, qh = 0.0, 0.0
    for w in target.wells
        pos = get_well_position(facility.domain, w)
        qw = state_facility.TotalSurfaceMassRate[pos]
        Tw = state_facility.SurfaceTemperature[pos]
        q += qw
        qh += qw.*Tw
    end
    T = qh./q

    return T
end

function well_top_node_enthalpy(ctrl::InjectorControl, model, state_well, T, cell)
    p = state_well.Pressure[cell]
    # density = ctrl.mixture_density
    # T = ctrl.temperature
    H_w = ctrl.enthalpy
    if ismissing(H_w)
        H = 0.0
        for ph in axes(state_well.Saturations, 1)
            # Define it via the volume weighted internal energy
            S = state_well.Saturations[ph, cell]
            dens = state_well.PhaseMassDensities[ph, cell]
            C = state_well.ComponentHeatCapacity[ph, cell]
            H += S*(C*T + p/dens)
        end
    elseif H_w isa Real
        H = H_w
    elseif H_w isa Function
        H = H_w(p, T)
    else
        error("InjectorControl.enthalpy must be missing, a real or a function (p, T).")
    end
    return H
end

function well_top_node_enthalpy(ctrl, model, state_well, T, cell)
    H = state_well.FluidEnthalpy
    S = state_well.Saturations
    H_w = 0.0
    for ph in axes(H, 1)
        H_w += H[ph, cell]*S[ph, cell]
    end
    return H_w
end

struct FacilityFromWellTemperatureCT <: Jutul.AdditiveCrossTerm
    well::Symbol
end

Jutul.cross_term_entities(ct::FacilityFromWellTemperatureCT, eq::SurfaceTemperatureEquation, model) = get_well_position(model.domain, ct.well)

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellTemperatureCT, eq, dt, ldisc = local_discretization(ct, i))

    pos = get_well_position(facility.domain, ct.well)
    T = 0*state_facility[:SurfaceTemperature][pos]
    T += state_well[:Temperature][well_top_node()]
    out[1] = -T
end