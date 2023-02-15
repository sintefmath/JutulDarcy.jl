export ReservoirFromWellFlowCT, FacilityFromWellFlowCT, WellFromFacilityFlowCT

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
    end
    rhoS = reference_densities(sys)

    p_well = state_s.Pressure
    p_res = state_t.Pressure
    # Wrap the key connection data in tuple for easy extension later
    conn = (dp = p_well[well_cell] - p_res[reservoir_cell],
           WI = WI, gdz = gdz,
           well = well_cell,
           reservoir = reservoir_cell)
    # Call smaller interface that is easy to specialize
    wg = model_s.domain.grid
    @inbounds well_perforation_flux!(out, wg, sys, state_t, state_s, rhoS, conn)
end

function perforation_phase_potential_difference(conn, state_res, state_well, ix)
    (; dp, WI, well, reservoir, gdz) = conn
    if gdz != 0
        ρ_r = state_res.PhaseMassDensities[ix, reservoir]
        if haskey(state_well, :PhaseMassDensities)
            ρ_w = state_well.PhaseMassDensities[ix, well]
            ρ = 0.5*(ρ_r + ρ_w)
        else
            ρ = ρ_r
        end
        dp += ρ*gdz
    end
    return -WI*dp
end

function Jutul.cross_term_entities(ct::AbstractReservoirFromWellCT, eq::Union{ConservationLaw, SimpleWellEquation}, model)
    return ct.reservoir_cells
end

function Jutul.cross_term_entities_source(ct::AbstractReservoirFromWellCT, eq::Union{ConservationLaw, SimpleWellEquation}, model)
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
            rhoS = reference_densities(well.system)
            rhoS, S = flash_wellstream_at_surface(well, state_well, rhoS)
            rhoS = tuple(rhoS...)
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

    if isa(ctrl, InjectorControl)
        if value(qT) < 0
            @warn "Injector $well_symbol is producing?"
        end
        mix = ctrl.injection_mixture
        nmix = length(mix)
        ncomp = number_of_components(flow_system(well.system))
        @assert nmix == ncomp "Injection composition length ($nmix) must match number of components ($ncomp)."
    else
        if value(qT) > 0
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
        WI = ct.WI[i]
    end

    p_well = state_well.Pressure[well_cell]
    p_res = state_res.Pressure[reservoir_cell]
    T_well = state_well.Temperature[well_cell]
    T_res = state_res.Temperature[reservoir_cell]

    kr = state_res.RelativePermeabilities
    mu = state_res.PhaseViscosities
    # Todo: Fix conn -> cell pressure drop
    ρgdz = 0
    mob = 0
    for ph in axes(kr, 1)
        mob += kr[ph, reservoir_cell]/mu[ph, reservoir_cell]
    end
    Q = -mob*WI*(p_well - p_res + ρgdz)
    fluid_heat = 0.0
    if Q < 0
        # Injection
        c = well_cell
        state_upw = state_well
    else
        c = reservoir_cell
        state_upw = state_res
    end
    S = state_upw.Saturations
    ρ = state_upw.PhaseMassDensities
    H = state_upw.FluidEnthalpy

    for ph in 1:nph
        fluid_heat += ρ[ph, c]*H[ph, c]*S[ph, c]
    end
    conductive_heat_flux = CI*(T_res - T_well)
    advective_heat_flux = fluid_heat*Q
    out[] = advective_heat_flux + conductive_heat_flux
end

function Base.show(io::IO, d::ReservoirFromWellThermalCT)
    n = length(d.CI)
    print(io, "ReservoirFromWellThermalCT ($n connections)")
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
    if isa(ctrl, InjectorControl)
        heat_capacity = state_well.FluidHeatCapacity[cell]
        p = state_well.Pressure[cell]
        density = ctrl.mixture_density
        T = ctrl.temperature
        H = heat_capacity*T + p/density
    else
        H = state_well.FluidEnthalpy[cell]
    end
    out[] = qT*H
end
