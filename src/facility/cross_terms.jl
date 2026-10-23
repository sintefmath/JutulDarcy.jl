
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
    apply_perforation_mask!(ct_s.target, mask, model.context)
    apply_perforation_mask!(ct_s.source, mask, model.context)
end

# Facility influence on well
struct WellFromFacilityFlowCT <: Jutul.AdditiveCrossTerm
    wells::Vector{Symbol}
    facility_cells::Vector{Int}
    well_cells::Vector{Int}
end

Jutul.cross_term_entities(ct::WellFromFacilityFlowCT, eq::ConservationLaw, model) =
    ct.well_cells

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
    well_symbol = ct.wells[i]
    pos = ct.facility_cells[i]
    cell = ct.well_cells[i]
    q_t, mix = cross_term_total_surface_mass_rate_and_mixture(
        facility, well, state_facility, state_well, well_symbol, pos, cell)
    for component in eachindex(out, mix)
        @inbounds out[component] = -mix[component]*q_t
    end
    return out
end

function cross_term_total_surface_mass_rate_and_mixture(facility, well, state_facility, state_well, well_symbol)
    pos = get_well_position(facility.domain, well_symbol)
    cell = WellMerging.well_top_node(physical_representation(well), well_symbol)
    return cross_term_total_surface_mass_rate_and_mixture(
        facility, well, state_facility, state_well, well_symbol, pos, cell)
end

function cross_term_total_surface_mass_rate_and_mixture(
        facility, well, state_facility, state_well, well_symbol, pos, cell)
    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    q_t = facility_surface_mass_rate_for_well(
        facility,
        well_symbol,
        state_facility,
        effective = true
    )
    # Hack for sparsity detection
    # bhp = bottom_hole_pressure(state_well)
    bhp = state_facility[:BottomHolePressure][pos]
    ph = state_facility[:TotalSurfaceMassRate][pos]
    ph2 = state_facility[:SurfacePhaseRates][1, pos]
    total_mass = state_well.TotalMasses[1, cell]
    q_t += 0*bhp + 0*ph + 0*total_mass + 0*ph2
    if isa(ctrl, InjectorControl)
        if value(q_t) < 0
            @warn "Injector $well_symbol is producing?"
        end
        mix = ctrl.injection_mixture
        nmix = length(mix)
        ncomp = number_of_components(well.system)
        @assert nmix == ncomp "Injection composition length ($nmix) must match number of components ($ncomp)."
    else
        if value(q_t) > 0 && ctrl isa ProducerControl
            @warn "Producer $well_symbol is injecting?"
        end
        if haskey(state_well, :MassFractions)
            mix = @view state_well.MassFractions[:, cell]
        else
            masses = @view state_well.TotalMasses[:, cell]
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
    if eq isa ConservationLaw{:TotalEnergy}
        # Account for the upwinded potential energy of the mass crossing the
        # perforation.
        qh += perforation_potential_energy_flux(λ_t, conn, state_res, state_well, nph)
    end
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

function perforation_potential_energy_flux(λ_t, conn, state_res, state_well, nph)
    q = 0
    for ph in 1:nph
        q_ph = perforation_phase_mass_flux(λ_t, conn, state_res, state_well, ph)
        if q_ph < 0
            # Injection
            Φ = state_well.UnitPotentialEnergy[conn.well]
        else
            Φ = state_res.UnitPotentialEnergy[conn.reservoir]
        end
        q += Φ*q_ph
    end
    return q
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
    apply_perforation_mask!(ct_s.target, mask, model.context)
    apply_perforation_mask!(ct_s.source, mask, model.context)
end

"Constitutive evaluator passed explicitly to the well/facility energy cross term."
struct InjectionEnthalpyEvaluator{E}
    energy::E
end

Adapt.@adapt_structure InjectionEnthalpyEvaluator

"""
Placeholder evaluator for fluid systems without a known way of computing the
injected enthalpy from pressure and temperature. Setup succeeds, so that
injectors with an explicitly given enthalpy can be used. An error is thrown only
if the injected enthalpy has to be computed.
"""
struct UnsupportedInjectionEnthalpy{E} end

"""
    injection_enthalpy_evaluator(model)

Return the evaluator used to compute the enthalpy of injected fluid from
pressure and temperature for injectors without an explicitly given enthalpy.
Packages with other fluid systems can add methods for their models, together
with a corresponding `injected_fluid_enthalpy` method.
"""
function injection_enthalpy_evaluator(model)
    # Resolve constitutive laws during simulator setup. Kernel models
    # intentionally do not carry variable definitions.
    enthalpy = model[:FluidEnthalpy]
    if enthalpy isa PressureTemperatureDependentEnthalpy
        return enthalpy
    end
    energy = model[:FluidInternalEnergy]
    if energy isa FluidInternalEnergy
        energy = nothing
    elseif !(energy isa TemperatureDependentVariable)
        return UnsupportedInjectionEnthalpy{typeof(energy)}()
    end
    return InjectionEnthalpyEvaluator(energy)
end

struct WellFromFacilityThermalCT{H} <: Jutul.AdditiveCrossTerm
    wells::Vector{Symbol}
    facility_cells::Vector{Int}
    well_cells::Vector{Int}
    enthalpy_evaluator::H
end

WellFromFacilityThermalCT(wells, facility_cells, well_cells) =
    WellFromFacilityThermalCT(wells, facility_cells, well_cells, nothing)

WellFromFacilityThermalCT(ct::WellFromFacilityThermalCT, evaluator) =
    WellFromFacilityThermalCT(ct.wells, ct.facility_cells, ct.well_cells, evaluator)

function Jutul.setup_cross_term_storage(ct::WellFromFacilityThermalCT{Nothing},
        eq_t, eq_s, well, facility, storage_t, storage_s; kwarg...)
    evaluator = injection_enthalpy_evaluator(well)
    storage = Jutul.setup_cross_term_storage(WellFromFacilityThermalCT(ct, evaluator),
        eq_t, eq_s, well, facility, storage_t, storage_s; kwarg...)
    storage[:injection_enthalpy_evaluator] = evaluator
    return storage
end

function Jutul.update_cross_term!(storage::Union{Jutul.JutulStorage, Jutul.ImmutableJutulStorage},
        ct::WellFromFacilityThermalCT{Nothing}, eq, storage_t, storage_s, well, facility, dt)
    cross_term = WellFromFacilityThermalCT(ct, storage.injection_enthalpy_evaluator)
    return Jutul.update_cross_term!(storage, cross_term, eq, storage_t, storage_s, well, facility, dt)
end

function Jutul.KernelExecution.setup_cross_term_evaluation(storage,
        ct::WellFromFacilityThermalCT{Nothing}, eq, storage_t, storage_s, well, facility)
    cross_term = WellFromFacilityThermalCT(ct, storage.injection_enthalpy_evaluator)
    return Jutul.KernelExecution.setup_cross_term_evaluation(storage, cross_term,
        eq, storage_t, storage_s, well, facility)
end

Jutul.cross_term_entities(ct::WellFromFacilityThermalCT, eq::ConservationLaw, model) =
    ct.well_cells

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
    well_symbol = ct.wells[i]
    pos = ct.facility_cells[i]

    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    qT = state_facility.TotalSurfaceMassRate[pos]
    # Hack for sparsity detection
    cell = ct.well_cells[i]
    qT += 0*state_well.Pressure[cell]

    H = get_target_enthalpy(ctrl, ctrl.target, facility, state_facility, well, state_well, cell, ct.enthalpy_evaluator)
    Q = -qT*H
    if eq isa ConservationLaw{:TotalEnergy}
        # Mass entering/leaving through the top node carries the potential
        # energy of that node
        Q -= qT*state_well.UnitPotentialEnergy[cell]
    end
    out[] = Q
end

function get_target_enthalpy(ctrl, target, facility, state_facility, model, state_well, cell, evaluator)
    return well_top_node_enthalpy(model, state_well, cell)
end

function get_target_enthalpy(ctrl::InjectorControl, target, facility, state_facility, model, state_well, cell, evaluator)
    T = get_target_temperature(ctrl, target, facility, state_facility)
    return well_top_node_enthalpy(ctrl, model, state_well, T, cell, evaluator)
end

function get_target_enthalpy(ctrl::InjectorControl, target::ReinjectionTarget, facility, state_facility, model, state_well, cell, evaluator)
    if !ismissing(ctrl.enthalpy) || !isnan(ctrl.temperature)
        T = get_target_temperature(ctrl, target, facility, state_facility)
        return well_top_node_enthalpy(ctrl, model, state_well, T, cell, evaluator)
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
        return state_well[:Enthalpy][cell]
    end
    H = state_well[:FluidEnthalpy]
    if haskey(state_well, :Saturations)
        S = state_well[:Saturations]
        rho = state_well[:PhaseMassDensities]
        H_w = zero(H[1, cell])
        total_mass = zero(H_w)
        for ph in axes(H, 1)
            mass = S[ph, cell]*rho[ph, cell]
            H_w += H[ph, cell]*mass
            total_mass += mass
        end
        H_w /= total_mass
    else
        H_w = H[1, cell]
    end
    return H_w
end

function well_top_node_enthalpy(ctrl::InjectorControl, model, state_well, T, cell, evaluator)
    p = state_well.Pressure[cell]
    H_w = ctrl.enthalpy
    if ismissing(H_w)
        H = injected_fluid_enthalpy(ctrl, evaluator, model, state_well, p, T, cell)
    elseif H_w isa Real
        H = H_w
    elseif H_w isa Function
        H = H_w(p, T)
    else
        error("InjectorControl.enthalpy must be missing, a real or a function (p, T).")
    end
    return H
end

function injected_fluid_enthalpy(ctrl, evaluator::InjectionEnthalpyEvaluator, model, state_well, p, T, cell)
    H = zero(p*T)
    law = evaluator.energy
    energy = evaluate_table_by_region(law.tab, region(law.regions, cell), T)
    for c in eachindex(ctrl.injection_mixture)
        H += ctrl.injection_mixture[c]*energy[c]
    end
    H += injection_pressure_work(ctrl, model, state_well, p, cell)
    return H
end

function injected_fluid_enthalpy(ctrl, evaluator::InjectionEnthalpyEvaluator{Nothing}, model, state, p, T, cell)
    H = zero(p*T)
    for c in eachindex(ctrl.injection_mixture)
        H += ctrl.injection_mixture[c]*state.ComponentHeatCapacity[c, cell]*T
    end
    return H + injection_pressure_work(ctrl, model, state, p, cell)
end

function injected_fluid_enthalpy(ctrl, ::UnsupportedInjectionEnthalpy{E}, model, state, p, T, cell) where E
    throw(ArgumentError("Cannot compute the injected enthalpy for $E. Give the \
        enthalpy explicitly in InjectorControl, or add methods for \
        injection_enthalpy_evaluator and injected_fluid_enthalpy."))
end

function injected_fluid_enthalpy(ctrl, enthalpy::PressureTemperatureDependentEnthalpy, model, state, p, T, cell)
    component_H = tabulated_enthalpy(enthalpy, p, T, cell)
    H = zero(p*T)
    for c in eachindex(ctrl.injection_mixture)
        H += ctrl.injection_mixture[c]*component_H[c]
    end
    return H
end

@inline function injection_pressure_work(ctrl, model, state, p, cell)
    work = zero(p)
    if model.system isa CompositionalSystem
        total_mass = zero(p)
        total_volume = zero(p)
        for ph in axes(state.Saturations, 1)
            saturation = state.Saturations[ph, cell]
            total_mass += saturation*state.PhaseMassDensities[ph, cell]
            total_volume += saturation
        end
        work = p*total_volume/total_mass
    else
        for ph in eachindex(ctrl.injection_mixture)
            work += ctrl.injection_mixture[ph]*p/state.PhaseMassDensities[ph, cell]
        end
    end
    return work
end

struct FacilityFromWellTemperatureCT <: Jutul.AdditiveCrossTerm
    wells::Vector{Symbol}
    facility_cells::Vector{Int}
    well_cells::Vector{Int}
end

Jutul.cross_term_entities(ct::FacilityFromWellTemperatureCT, eq::SurfaceTemperatureEquation, model) = ct.facility_cells

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellTemperatureCT, eq, dt, ldisc = local_discretization(ct, i))

    pos = ct.facility_cells[i]
    cell = ct.well_cells[i]
    T = 0*state_facility[:SurfaceTemperature][pos]
    T += state_well[:Temperature][cell]
    out[1] = -T
end

struct FacilityFromWellEnthalpyCT <: Jutul.AdditiveCrossTerm
    wells::Vector{Symbol}
    facility_cells::Vector{Int}
    well_cells::Vector{Int}
end

Jutul.cross_term_entities(ct::FacilityFromWellEnthalpyCT, eq::SurfaceEnthalpyEquation, model) = ct.facility_cells

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellEnthalpyCT, eq::SurfaceEnthalpyEquation, dt, ldisc = local_discretization(ct, i))

    pos = ct.facility_cells[i]
    cell = ct.well_cells[i]
    H = 0*state_facility[:SurfaceEnthalpy][pos]
    H += well_top_node_enthalpy(well, state_well, cell)
    out[1] = -H*eq.scale
end

struct FacilityFromWellBottomHolePressureCT <: Jutul.AdditiveCrossTerm
    wells::Vector{Symbol}
    facility_cells::Vector{Int}
    well_cells::Vector{Int}
end

Jutul.cross_term_entities(ct::FacilityFromWellBottomHolePressureCT, eq::BottomHolePressureEquation, model) = ct.facility_cells

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellBottomHolePressureCT, eq::BottomHolePressureEquation, dt, ldisc = local_discretization(ct, i))

    pos = ct.facility_cells[i]
    cell = ct.well_cells[i]
    P = 0*state_facility[:BottomHolePressure][pos]
    P += state_well[:Pressure][cell]
    out[1] = -P*eq.scale
end

struct FacilityFromWellSurfaceComponentRatesCT <: Jutul.AdditiveCrossTerm
    wells::Vector{Symbol}
    facility_cells::Vector{Int}
    well_cells::Vector{Int}
end

Jutul.cross_term_entities(ct::FacilityFromWellSurfaceComponentRatesCT, eq::SurfaceComponentRatesEquation, model) = ct.facility_cells

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    ct::FacilityFromWellSurfaceComponentRatesCT, eq::SurfaceComponentRatesEquation, dt, ldisc = local_discretization(ct, i))

    well_symbol = ct.wells[i]
    pos = ct.facility_cells[i]
    cell = ct.well_cells[i]
    q_t, mixture = cross_term_total_surface_mass_rate_and_mixture(
        facility, well, state_facility, state_well, well_symbol, pos, cell)
    for component in eachindex(out, mixture)
        out[component] = -q_t*mixture[component]*eq.scale
    end
    return out
end
