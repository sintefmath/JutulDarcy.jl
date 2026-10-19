struct WellFromFacilityTracerCT{T} <: Jutul.AdditiveCrossTerm
    wells::Vector{Symbol}
    facility_cells::Vector{Int}
    well_cells::Vector{Int}
    tracers::T
end

Jutul.cross_term_entities(ct::WellFromFacilityTracerCT, eq::ConservationLaw{:TracerMasses}, model) =
    ct.well_cells

const TRACER_TOL = 1e-12

function Jutul.update_cross_term_in_entity!(out, i,
    state_well, state0_well,
    state_facility, state0_facility,
    well, facility,
    ct::WellFromFacilityTracerCT, eq, dt, ldisc = Jutul.local_discretization(ct, i))
    well_symbol = ct.wells[i]
    pos = ct.facility_cells[i]

    cfg = state_facility.WellGroupConfiguration
    ctrl = JutulDarcy.operating_control(cfg, well_symbol)
    qT = state_facility.TotalSurfaceMassRate[pos]

    tracers = ct.tracers
    T = eltype(out)
    N = length(tracers)
    S = state_well.Saturations
    rho = state_well.PhaseMassDensities
    wc = ct.well_cells[i]

    if ctrl isa InjectorControl
        if ismissing(ctrl.tracers)
            out .= zero(T)
        else
            @assert length(ctrl.tracers) == N
            update_tracer_injection!(out, ctrl, qT, tracers, well.system, state_well, wc)
        end
    else
        mass_tot = zero(T)
        for phase in 1:number_of_phases(well.system)
            mass_tot += S[phase, wc]*rho[phase, wc]
        end
        for i in 1:N
            tracer = tracers[i]
            v = zero(T)
            C_i = state_well.TracerConcentrations[i, wc]
            for phase in tracer_phase_indices(tracer)
                F_ph = S[phase, wc]*rho[phase, wc]/mass_tot
                v += F_ph*C_i*qT
            end
            out[i] = -v
        end
    end
    return out
end

function update_tracer_injection!(out, ctrl, qT, tracers,
        system::Union{JutulDarcy.BlackOilSystem, ImmiscibleSystem}, state, cell)
    for i in eachindex(tracers)
        F = zero(eltype(out))
        for phase in tracer_phase_indices(tracers[i])
            F += ctrl.injection_mixture[phase]
        end
        out[i] = -qT*F*ctrl.tracers[i]
    end
    return out
end

function update_tracer_injection!(out, ctrl, qT, tracers, system, state, cell)
    S = state.Saturations
    rho = state.PhaseMassDensities
    mass_tot = zero(eltype(out))
    for phase in 1:number_of_phases(system)
        mass_tot += S[phase, cell]*rho[phase, cell]
    end
    for i in eachindex(tracers)
        mass = zero(eltype(out))
        for phase in tracer_phase_indices(tracers[i])
            mass += S[phase, cell]*rho[phase, cell]
        end
        if mass > TRACER_TOL
            out[i] = -qT*(mass/mass_tot)*ctrl.tracers[i]
        else
            out[i] = zero(eltype(out))
        end
    end
    return out
end

struct ReservoirFromWellTracerCT{I<:AbstractVector} <: JutulDarcy.AbstractReservoirFromWellCT
    reservoir_cells::I
    well_cells::I
end

function Adapt.adapt_structure(to, ct::ReservoirFromWellTracerCT)
    return ReservoirFromWellTracerCT(
        Adapt.adapt(to, ct.reservoir_cells),
        Adapt.adapt(to, ct.well_cells))
end

function Jutul.update_cross_term_in_entity!(out, i,
    state_res, state0_res,
    state_well, state0_well, 
    model_res, model_well,
    ct::ReservoirFromWellTracerCT, eq, dt, ldisc = Jutul.local_discretization(ct, i))
    # Unpack properties
    sys = model_res.system
    nph = number_of_phases(sys)
    @inbounds begin 
        reservoir_cell = ct.reservoir_cells[i]
        well_cell = ct.well_cells[i]
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
    λ_t = zero(eltype(state_res.PhaseMobilities))
    for ph in 1:nph
        λ_t += state_res.PhaseMobilities[ph, reservoir_cell]
    end

    tracers = eq.flux_type.tracers
    T = eltype(out)
    for tracer_index in eachindex(tracers)
        tracer = tracers[tracer_index]
        v = zero(T)
        for phase in tracer_phase_indices(tracer)
            q_ph = JutulDarcy.perforation_phase_mass_flux(
                λ_t, conn, state_res, state_well, phase)
            if q_ph < 0
                # Injection
                C_perf = state_well.TracerConcentrations[
                    tracer_index, well_cell]
            else
                C_perf = state_res.TracerConcentrations[
                    tracer_index, reservoir_cell]
            end
            v += C_perf*q_ph
        end
        out[tracer_index] = v
    end
    return out
end

function Base.show(io::IO, d::ReservoirFromWellTracerCT)
    n = length(d.reservoir_cells)
    print(io, "ReservoirFromWellTracerCT ($n connections)")
end

function Jutul.subcrossterm(ct::ReservoirFromWellTracerCT, ctp, m_t, m_s, map_res::Jutul.FiniteVolumeGlobalMap, ::Jutul.TrivialGlobalMap, partition)
    (; reservoir_cells, well_cells) = ct
    rc = map(
        c -> Jutul.local_cell(c, map_res),
        reservoir_cells)
    return ReservoirFromWellTracerCT(rc, copy(well_cells))
end

function Jutul.apply_force_to_cross_term!(ct_s, cross_term::ReservoirFromWellTracerCT, target, source, model, storage, dt, force::PerforationMask; time = time)
    mask = force.values
    JutulDarcy.apply_perforation_mask!(ct_s.target, mask, model.context)
    JutulDarcy.apply_perforation_mask!(ct_s.source, mask, model.context)
end
