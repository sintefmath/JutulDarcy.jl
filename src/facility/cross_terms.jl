
abstract type AbstractReservoirFromWellCT <: Jutul.AdditiveCrossTerm end
struct ReservoirFromWellFlowCT{I<:AbstractVector} <: AbstractReservoirFromWellCT
    reservoir_cells::I
    well_cells::I
end

function Base.show(io::IO, d::ReservoirFromWellFlowCT)
    n = length(d.reservoir_cells)
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
    return out
end

facility_raw_state_field(state, key::Symbol) = state[key]
facility_raw_state_field(state::Union{Jutul.LocalStateAD, Jutul.ValueStateAD}, key::Symbol) =
    getfield(state, :data)[key]

facility_is_injector(::Nothing, state, well, pos) =
    @inbounds facility_raw_state_field(state, :FacilityCrossTermState).control_type[pos] == 1
facility_is_injector(cfg, state, well, pos) =
    operating_control(cfg, well) isa InjectorControl

facility_control_factor(::Nothing, state, well, pos) =
    @inbounds facility_raw_state_field(state, :FacilityCrossTermState).factor[pos]
function facility_control_factor(cfg, state, well, pos)
    control = operating_control(cfg, well)
    return control isa Union{InjectorControl, ProducerControl} ? control.factor : 1.0
end

facility_injection_component(::Nothing, state, well, pos, component) =
    @inbounds facility_raw_state_field(state, :FacilityCrossTermState).injection_mixture[component, pos]
facility_injection_component(cfg, state, well, pos, component) =
    @inbounds operating_control(cfg, well).injection_mixture[component]

facility_injection_density(::Nothing, state, well, pos) =
    @inbounds facility_raw_state_field(state, :FacilityCrossTermState).mixture_density[pos]
facility_injection_density(cfg, state, well, pos) =
    operating_control(cfg, well).mixture_density

facility_injection_phase_fraction(::Nothing, state, well, pos, phase) =
    @inbounds facility_raw_state_field(state, :FacilityCrossTermState).phase_fractions[phase, pos]
function facility_injection_phase_fraction(cfg, state, well, pos, phase)
    fraction = 0.0
    for (phase_index, phase_fraction) in operating_control(cfg, well).phases
        if phase_index == phase
            fraction = phase_fraction
        end
    end
    return fraction
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

Base.@propagate_inbounds function perforation_phase_potential_difference(conn, state_res, state_well, ix::Int)
    dp = conn.dp
    WI = conn.WI
    WI, dp = Base.promote(WI, dp)
    if haskey(state_res, :PermeabilityMultiplier)
        K_mul = state_res[:PermeabilityMultiplier][conn.reservoir]
        WI *= K_mul
    end
    if conn.gdz != 0.0
        if haskey(state_well, :ConnectionPressureDrop)
            cdp = state_well.ConnectionPressureDrop[conn.perforation]
            dp += cdp
        else
            ρ_r = state_res.PhaseMassDensities[ix, conn.reservoir]
            if haskey(state_well, :PhaseMassDensities)
                ρ_w = state_well.PhaseMassDensities[ix, conn.well]
                ρ = 0.5*(ρ_r + ρ_w)
            else
                ρ = ρ_r
            end
            dp += ρ*conn.gdz
        end
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
    (; reservoir_cells, well_cells) = ct
    rc = map(
        c -> Jutul.local_cell(c, map_res),
        reservoir_cells)
    return ReservoirFromWellFlowCT(rc, copy(well_cells))
end

function well_top_node()
    return 1
end

function Jutul.apply_force_to_cross_term!(ct_s, cross_term::ReservoirFromWellFlowCT, target, source, model, storage, dt, force::PerforationMask; time = time)
    mask = force.values
    apply_perforation_mask!(ct_s.target, mask)
    apply_perforation_mask!(ct_s.source, mask)
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
struct WellFromFacilityFlowCT{W} <: Jutul.AdditiveCrossTerm
    well::W
    facility_position::Int
end
WellFromFacilityFlowCT(well) = WellFromFacilityFlowCT(well, 0)

@inline function facility_position(ct, facility)
    position = ct.facility_position
    return iszero(position) ? get_well_position(facility.domain, ct.well) : position
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
        ct::WellFromFacilityFlowCT, eq, dt, ldisc = local_discretization(ct, i)
    )
    well_symbol = ct.well
    cross_term_surface_mass_rate!(
        out, facility, well, state_facility, state_well, well_symbol,
        facility_position(ct, facility))
    return out
end

function cross_term_surface_mass_rate!(out, facility, well, state_facility,
        state_well, well_symbol, pos = get_well_position(
            facility.domain, well_symbol))
    cfg = facility_raw_state_field(state_facility, :WellGroupConfiguration)
    q_t = state_facility.TotalSurfaceMassRate[pos]
    q_t *= facility_control_factor(cfg, state_facility, well_symbol, pos)

    # Preserve all intended AD dependencies even when the active control makes
    # some of these values numerically inactive.
    bhp = state_facility.BottomHolePressure[pos]
    phase_rate = state_facility.SurfacePhaseRates[1, pos]
    total_mass = state_well.TotalMasses[1, well_top_node()]
    q_t += 0*(bhp + phase_rate + total_mass)

    if facility_is_injector(cfg, state_facility, well_symbol, pos)
        for component in eachindex(out)
            fraction = facility_injection_component(
                cfg, state_facility, well_symbol, pos, component)
            out[component] = -fraction*q_t
        end
    elseif haskey(state_well, :MassFractions)
        fractions = state_well.MassFractions
        for component in eachindex(out)
            out[component] = -fractions[component, well_top_node()]*q_t
        end
    else
        masses = state_well.TotalMasses
        mass = zero(masses[1, well_top_node()])
        for component in eachindex(out)
            mass += masses[component, well_top_node()]
        end
        for component in eachindex(out)
            out[component] = -(masses[component, well_top_node()]/mass)*q_t
        end
    end
    return out
end

function cross_term_total_surface_mass_rate_and_mixture(facility, well, state_facility, state_well, well_symbol)
    pos = get_well_position(facility.domain, well_symbol)

    cfg = facility_raw_state_field(state_facility, :WellGroupConfiguration)
    q_t = state_facility.TotalSurfaceMassRate[pos]
    q_t *= facility_control_factor(cfg, state_facility, well_symbol, pos)
    # Hack for sparsity detection
    # bhp = bottom_hole_pressure(state_well)
    bhp = state_facility[:BottomHolePressure][pos]
    ph = state_facility[:TotalSurfaceMassRate][pos]
    ph2 = state_facility[:SurfacePhaseRates][1, pos]
    total_mass = state_well.TotalMasses[1, well_top_node()]
    q_t += 0*bhp + 0*ph + 0*total_mass + 0*ph2
    if facility_is_injector(cfg, state_facility, well_symbol, pos)
        nmix = size(facility_raw_state_field(
            state_facility, :FacilityCrossTermState).injection_mixture, 1)
        ncomp = number_of_components(well.system)
        @assert nmix == ncomp "Injection composition length ($nmix) must match number of components ($ncomp)."
        mix = ntuple(component -> facility_injection_component(
            cfg, state_facility, well_symbol, pos, component), ncomp)
    else
        if haskey(state_well, :MassFractions)
            mix = state_well.MassFractions
        else
            masses = @views state_well.TotalMasses[:, well_top_node()]
            mass = sum(masses)
            mix = masses./mass
        end
    end
    return (q_t, mix)
end

# Thermal
struct ReservoirFromWellThermalCT{I<:AbstractVector} <: AbstractReservoirFromWellCT
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
        WI = state_well.WellIndices[i]
        WIth = state_well.WellIndicesThermal[i]
        gdz = state_well.PerforationGravityDifference[i]
        p_well = state_well.Pressure
        p_res = state_res.Pressure
        dp = p_well[well_cell] - p_res[reservoir_cell]
        conn = (
            dp = dp,
            WI = WI,
            WIth = WIth,
            gdz = gdz,
            well = well_cell,
            perforation = i,
            reservoir = reservoir_cell
        )
    end

    λ_t = sum(perforation_reservoir_mobilities(state_res, state_well, sys, reservoir_cell, well_cell))
    qh = perforation_phase_thermal_flux(λ_t, conn, state_res, state_well, nph)
    out[] = qh

end

function perforation_phase_thermal_flux(λ_t, conn, state_res, state_well, nph)

    well_cell = conn.well
    reservoir_cell = conn.reservoir
    WIth = conn.WIth

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
    return advective_heat_flux + conductive_heat_flux
end

function Base.show(io::IO, d::ReservoirFromWellThermalCT)
    n = length(d.well_cells)
    print(io, "ReservoirFromWellThermalCT ($n connections)")
end

function Jutul.subcrossterm(ct::ReservoirFromWellThermalCT, ctp, m_t, m_s, map_res::FiniteVolumeGlobalMap, ::TrivialGlobalMap, partition)
    (; reservoir_cells, well_cells) = ct
    rc = map(
        c -> Jutul.local_cell(c, map_res),
        reservoir_cells)
    return ReservoirFromWellThermalCT(rc, copy(well_cells))
end

function Jutul.apply_force_to_cross_term!(ct_s, cross_term::ReservoirFromWellThermalCT, target, source, model, storage, dt, force::PerforationMask; time = time)
    mask = force.values
    apply_perforation_mask!(ct_s.target, mask)
    apply_perforation_mask!(ct_s.source, mask)
end

struct WellFromFacilityThermalCT{W} <: Jutul.AdditiveCrossTerm
    well::W
    facility_position::Int
end
WellFromFacilityThermalCT(well) = WellFromFacilityThermalCT(well, 0)

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
    pos = facility_position(ct, facility)

    cfg = facility_raw_state_field(state_facility, :WellGroupConfiguration)
    ctrl = operating_control(cfg, well_symbol)
    qT = state_facility.TotalSurfaceMassRate[pos]
    # Hack for sparsity detection
    qT += 0*bottom_hole_pressure(state_well)

    cell = well_top_node()

    H = get_target_enthalpy(ctrl, ctrl.target, facility, state_facility, well, state_well, cell)
    out[] = -qT*H
end

function get_target_enthalpy(ctrl, target, facility, state_facility, model, state_well, cell)
    return well_top_node_enthalpy(model, state_well, cell)
end

function get_target_enthalpy(ctrl::InjectorControl, target, facility, state_facility, model, state_well, cell)
    T = get_target_temperature(ctrl, target, facility, state_facility)
    return well_top_node_enthalpy(ctrl, model, state_well, T, cell)
end

function get_target_enthalpy(ctrl::InjectorControl, target::ReinjectionTarget, facility, state_facility, model, state_well, cell)
    if !ismissing(ctrl.enthalpy) || !isnan(ctrl.temperature)
        T = get_target_temperature(ctrl, target, facility, state_facility)
        return well_top_node_enthalpy(ctrl, model, state_well, T, cell)
    end

    q = qh = Htot = 0.0
    for w in target.wells
        pos = get_well_position(facility.domain, w)
        qw = state_facility.TotalSurfaceMassRate[pos]
        Hw = state_facility.SurfaceEnthalpy[pos]
        q += qw
        qh += qw*Hw
        Htot += Hw
    end
    return ifelse(abs(q) >= MIN_ACTIVE_WELL_RATE, qh/q, Htot/length(target.wells))
end

function get_target_temperature(ctrl, target, facility, state_facility)
    return missing
end

function get_target_temperature(ctrl::InjectorControl, target, facility, state_facility)
    return ctrl.temperature
end

function get_target_temperature(ctrl::InjectorControl, target::ReinjectionTarget, facility, state_facility)
    if !isnan(ctrl.temperature)
        return ctrl.temperature
    end

    q = qh = Ttot = 0.0
    for w in target.wells
        pos = get_well_position(facility.domain, w)
        qw = state_facility.TotalSurfaceMassRate[pos]
        Tw = state_facility.SurfaceTemperature[pos]
        q += qw
        qh += qw.*Tw
        Ttot += Tw
    end
    T = ifelse(abs(q) >= MIN_ACTIVE_WELL_RATE, qh/q, Ttot/length(target.wells))

    return T
end

function well_top_node_enthalpy(model, state_well, cell)
    if haskey(state_well, :Enthalpy)
        return state_well.Enthalpy[cell]
    end
    H = state_well.FluidEnthalpy
    if haskey(state_well, :Saturations)
        S = state_well.Saturations
        H_w = zero(H[1, cell])
        for ph in axes(H, 1)
            H_w += H[ph, cell]*S[ph, cell]
        end
    else
        H_w = H[1, cell]
    end
    return H_w
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
    return well_top_node_enthalpy(model, state_well, cell)
end

struct FacilityFromWellTemperatureCT{W} <: Jutul.AdditiveCrossTerm
    well::W
    facility_position::Int
end
FacilityFromWellTemperatureCT(well) = FacilityFromWellTemperatureCT(well, 0)

Jutul.cross_term_entities(ct::FacilityFromWellTemperatureCT, eq::SurfaceTemperatureEquation, model) = facility_position(ct, model)

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellTemperatureCT, eq, dt, ldisc = local_discretization(ct, i))

    pos = facility_position(ct, facility)
    T = 0*state_facility[:SurfaceTemperature][pos]
    T += state_well[:Temperature][well_top_node()]
    out[1] = -T
end

struct FacilityFromWellEnthalpyCT{W} <: Jutul.AdditiveCrossTerm
    well::W
    facility_position::Int
end
FacilityFromWellEnthalpyCT(well) = FacilityFromWellEnthalpyCT(well, 0)

Jutul.cross_term_entities(ct::FacilityFromWellEnthalpyCT, eq::SurfaceEnthalpyEquation, model) = facility_position(ct, model)

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellEnthalpyCT, eq::SurfaceEnthalpyEquation, dt, ldisc = local_discretization(ct, i))

    pos = facility_position(ct, facility)
    H = 0*state_facility[:SurfaceEnthalpy][pos]
    H += well_top_node_enthalpy(well, state_well, well_top_node())
    out[1] = -H*eq.scale
end

struct FacilityFromWellBottomHolePressureCT{W} <: Jutul.AdditiveCrossTerm
    well::W
    facility_position::Int
end
FacilityFromWellBottomHolePressureCT(well) =
    FacilityFromWellBottomHolePressureCT(well, 0)

Jutul.cross_term_entities(ct::FacilityFromWellBottomHolePressureCT, eq::BottomHolePressureEquation, model) = facility_position(ct, model)

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellBottomHolePressureCT, eq::BottomHolePressureEquation, dt, ldisc = local_discretization(ct, i))

    pos = facility_position(ct, facility)
    P = 0*state_facility[:BottomHolePressure][pos]
    P += state_well[:Pressure][well_top_node()]
    out[1] = -P*eq.scale
end

struct FacilityFromSurfacePhaseRatesCT{W} <: Jutul.AdditiveCrossTerm
    well::W
    facility_position::Int
end
FacilityFromSurfacePhaseRatesCT(well) = FacilityFromSurfacePhaseRatesCT(well, 0)

const DeviceFacilityCrossTerm = Union{
    WellFromFacilityFlowCT{Nothing},
    WellFromFacilityThermalCT{Nothing},
    FacilityFromWellTemperatureCT{Nothing},
    FacilityFromWellEnthalpyCT{Nothing},
    FacilityFromWellBottomHolePressureCT{Nothing},
    FacilityFromSurfacePhaseRatesCT{Nothing}
}

# Adapted device cross terms have no well name and always carry the resolved
# integer position. This dispatch keeps that position concrete in GPU kernels.
@inline facility_position(ct::DeviceFacilityCrossTerm, facility) =
    ct.facility_position

Jutul.cross_term_entities(ct::FacilityFromSurfacePhaseRatesCT, eq::SurfacePhaseRatesEquation, model) = facility_position(ct, model)

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromSurfacePhaseRatesCT, eq::SurfacePhaseRatesEquation, dt, ldisc = local_discretization(ct, i))

    pos = facility_position(ct, facility)
    q_t = state_facility.TotalSurfaceMassRate[pos]
    cfg = facility_raw_state_field(state_facility, :WellGroupConfiguration)
    rhoS, S = surface_density_and_volume_fractions(state_well)

    p_top = state_well.Pressure[well_top_node()]
    tm = state_well.TotalMasses[1, well_top_node()]
    # Sparsity hack
    for ph_idx in eachindex(out)
        out[ph_idx] = 0.0*(S[ph_idx] + rhoS[ph_idx] + q_t + tm + p_top)
    end
    is_injector = facility_is_injector(cfg, state_facility, ct.well, pos)
    if is_injector
        density = facility_injection_density(cfg, state_facility, ct.well, pos)
        volume_rate = q_t/density
        for ph_idx in eachindex(out)
            fraction = facility_injection_phase_fraction(
                cfg, state_facility, ct.well, pos, ph_idx)
            out[ph_idx] += -volume_rate*fraction*eq.scale
        end
    else
        total_density = 0.0
        for i in eachindex(rhoS, S)
            total_density += S[i]*rhoS[i]
        end
        q_vol = q_t/total_density
        for ph in eachindex(rhoS, S)
            out[ph] += -S[ph]*q_vol*eq.scale
        end
    end
    return out
end
