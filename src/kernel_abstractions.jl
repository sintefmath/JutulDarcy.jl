# Structural adaptation for the compact topology used by the single-model
# reservoir discretization. More complex property types can add equivalent
# Adapt rules without changing the execution context.
function Adapt.adapt_structure(to, g::MinimalTPFATopology)
    return MinimalTPFATopology(g.nc, Adapt.adapt(to, g.neighborship))
end

Adapt.@adapt_structure Rs
Adapt.@adapt_structure Rv
Adapt.@adapt_structure PVTO
Adapt.@adapt_structure PVTOTable
Adapt.@adapt_structure PVDO
Adapt.@adapt_structure PVDG
Adapt.@adapt_structure PVTG
Adapt.@adapt_structure PVTGTable
Adapt.@adapt_structure PVTW
Adapt.@adapt_structure FacilitySystem

function Adapt.adapt_structure(to, table::DeckThermalViscosityTable)
    # The outer vector is a small set of PVT regions. Keeping it as a tuple
    # lets Adapt recursively transfer each interpolant's arrays.
    visc_tab = map(table.visc_tab) do phase_tables
        map(x -> Adapt.adapt(to, x), Tuple(phase_tables))
    end
    return DeckThermalViscosityTable(
        visc_tab,
        Adapt.adapt(to, table.p_ref),
        Adapt.adapt(to, table.rs_ref))
end

function Adapt.adapt_structure(to,
        variable::TemperatureDependentVariable{T, R, N}) where {T, R, N}
    return TemperatureDependentVariable(
        Adapt.adapt(to, variable.tab),
        Adapt.adapt(to, variable.regions),
        Val(N))
end

function Adapt.adapt_structure(to, variable::SimpleCapillaryPressure)
    return SimpleCapillaryPressure(
        Adapt.adapt(to, variable.pc),
        Adapt.adapt(to, variable.regions),
        Val(:assembled))
end

function Adapt.adapt_structure(to, variable::ScaledCapillaryPressure)
    return ScaledCapillaryPressure(
        Adapt.adapt(to, variable.pc),
        Adapt.adapt(to, variable.regions),
        Val(:assembled))
end

function Adapt.adapt_structure(::Jutul.KernelAbstractionsContext, group::WellGroup)
    # The authoritative WellGroup, including its symbol vector, remains on the
    # host. Device assembly only needs the number of facility unknowns; cross
    # terms carry their precomputed integer facility positions.
    return WellGroup(
        Base.OneTo(length(group.well_symbols)),
        group.can_shut_producers,
        group.can_shut_injectors)
end

function Adapt.adapt_structure(to, well::SimpleWell)
    return SimpleWell(
        Adapt.adapt(to, well.perforations),
        Adapt.adapt(to, well.surface),
        nothing,
        well.explicit_dp
    )
end

function Adapt.adapt_structure(to,
        system::StandardBlackOilSystem{D, V, W, R, F}) where {D, V, W, R, F}
    rs_max = Adapt.adapt(to, system.rs_max)
    rv_max = Adapt.adapt(to, system.rv_max)
    rho_ref = Adapt.adapt(to, system.rho_ref)
    phase_indices = Adapt.adapt(to, system.phase_indices)
    phases = Adapt.adapt(to, system.phases)
    Num = typeof(system.rs_eps)
    return StandardBlackOilSystem{
        typeof(rs_max), typeof(rv_max), W, typeof(rho_ref), F,
        typeof(phase_indices), typeof(phases), Num
    }(rs_max, rv_max, rho_ref, phase_indices, phases,
        system.saturated_chop, system.keep_bubble_flag,
        system.rs_eps, system.rv_eps, system.s_eps,
        system.reference_phase_index)
end

function Adapt.adapt_structure(to,
        variable::AbstractReservoirRelativePermeabilities{Scaling, ph}) where {Scaling, ph}
    krw = Adapt.adapt(to, variable.krw)
    krow = Adapt.adapt(to, variable.krow)
    krog = Adapt.adapt(to, variable.krog)
    krg = Adapt.adapt(to, variable.krg)
    regions = Adapt.adapt(to, variable.regions)
    hysteresis_w = Adapt.adapt(to, variable.hysteresis_w)
    hysteresis_ow = Adapt.adapt(to, variable.hysteresis_ow)
    hysteresis_og = Adapt.adapt(to, variable.hysteresis_og)
    hysteresis_g = Adapt.adapt(to, variable.hysteresis_g)
    scaling = Adapt.adapt(to, variable.scaling)
    method = Adapt.adapt(to, variable.three_phase_method)
    return BackendReservoirRelativePermeabilities{
        typeof(scaling), ph, typeof(krw), typeof(krow), typeof(krog), typeof(krg),
        typeof(regions), typeof(hysteresis_w), typeof(hysteresis_ow),
        typeof(hysteresis_og), typeof(hysteresis_g), typeof(method)
    }(krw, krow, krog, krg, regions,
        hysteresis_w, hysteresis_ow, hysteresis_og, hysteresis_g, scaling,
        variable.hysteresis_s_threshold, variable.hysteresis_s_eps, method)
end

function Adapt.adapt_structure(to, variable::PhaseRelativePermeability)
    return BackendPhaseRelativePermeability(
        Adapt.adapt(to, variable.k), Val(variable.label),
        variable.connate, variable.critical, variable.s_max,
        variable.k_max, variable.input_s_max)
end

function Adapt.adapt_structure(to,
        variable::BackendPhaseRelativePermeability{label}) where label
    return BackendPhaseRelativePermeability(
        Adapt.adapt(to, variable.k), Val(label),
        variable.connate, variable.critical, variable.s_max,
        variable.k_max, variable.input_s_max)
end

function Adapt.adapt_structure(to, table::MuBTable)
    return MuBTable(
        Adapt.adapt(to, table.pressure),
        Adapt.adapt(to, table.shrinkage),
        Adapt.adapt(to, table.shrinkage_interp),
        Adapt.adapt(to, table.viscosity),
        Adapt.adapt(to, table.viscosity_interp),
        Val(:assembled)
    )
end

function Adapt.adapt_structure(to, variable::DeckPhaseViscosities)
    return DeckPhaseViscosities(
        Adapt.adapt(to, variable.pvt),
        Adapt.adapt(to, variable.thermal),
        Adapt.adapt(to, variable.regions),
        Val(:assembled))
end

function Adapt.adapt_structure(to, variable::DeckPhaseMassDensities)
    return DeckPhaseMassDensities(
        Adapt.adapt(to, variable.pvt),
        Adapt.adapt(to, variable.watdent),
        Adapt.adapt(to, variable.regions),
        Val(:assembled))
end

function Adapt.adapt_structure(to, variable::DeckShrinkageFactors)
    return DeckShrinkageFactors(
        Adapt.adapt(to, variable.pvt),
        Adapt.adapt(to, variable.watdent),
        Adapt.adapt(to, variable.regions),
        Val(:assembled))
end

function Adapt.adapt_structure(to, variable::LinearlyCompressiblePoreVolume)
    return LinearlyCompressiblePoreVolume(
        Adapt.adapt(to, variable.reference_pressure),
        Adapt.adapt(to, variable.expansion),
        Adapt.adapt(to, variable.regions),
        Val(:assembled))
end

struct BackendSurfaceWellConditions <: Jutul.ScalarVariable end

function Adapt.adapt_structure(::Jutul.KernelAbstractionsContext,
        variable::SurfaceWellConditions)
    isempty(variable.separator_conditions) || throw(ArgumentError(
        "GPU well execution does not support separator stages"))
    return BackendSurfaceWellConditions()
end

Jutul.get_dependencies(::BackendSurfaceWellConditions, model) = [:TotalMasses]

function Jutul.update_secondary_variable!(x::AbstractVector{TopConditions{N, R}},
        ::BackendSurfaceWellConditions, model, state, ix) where {N, R}
    rho = reference_densities(model.system)
    masses = haskey(state, :MassFractions) ? state.MassFractions : state.TotalMasses
    surface_volume = ntuple(Val(N)) do component
        @inbounds masses[component, 1]/rho[component]
    end
    total_volume = sum(surface_volume)
    fractions = ntuple(Val(N)) do component
        surface_volume[component]/total_volume
    end
    @inbounds x[1] = TopConditions(Val(N), Val(R), rho, fractions)
    return x
end

# The dictionary-based control configuration is host-only. Device cross terms
# consume FacilityCrossTermState, which is refreshed and copied in place after
# every host facility update.
Adapt.adapt_structure(::Jutul.KernelAbstractionsContext,
    ::WellGroupConfiguration) = nothing

function Adapt.adapt_structure(to, state::FacilityCrossTermState)
    return FacilityCrossTermState(
        Adapt.adapt(to, state.control_type),
        Adapt.adapt(to, state.factor),
        Adapt.adapt(to, state.mixture_density),
        Adapt.adapt(to, state.injection_temperature),
        Adapt.adapt(to, state.injection_enthalpy),
        Adapt.adapt(to, state.injection_mixture),
        Adapt.adapt(to, state.phase_fractions),
        Adapt.adapt(to, state.tracer_concentrations)
    )
end

function Jutul.backend_copyto!(destination::FacilityCrossTermState,
        source::FacilityCrossTermState)
    return copyto!(destination, source)
end

# Application extensions can report additional numeric facility data that must
# be allocated before the storage is adapted. Tracers uses this hook without
# making the general reservoir/well transfer depend on the Tracers module.
backend_facility_number_of_tracers(model) = 0

function Jutul.prepare_backend_transfer!(storage, model::Jutul.MultiModel)
    T = Jutul.float_type(model.context)
    for (index, pair) in enumerate(model.cross_terms)
        cross_term = pair.cross_term
        if cross_term isa AbstractReservoirFromWellCT
            cross_term_storage = storage.cross_terms[index]
            force_buffer = zeros(T, length(cross_term.reservoir_cells))
            storage.cross_terms[index] = Jutul.JutulStorage(merge(
                Jutul.data(cross_term_storage), (; force_buffer)))
        end
    end

    ntracers = backend_facility_number_of_tracers(model)
    iszero(ntracers) && return storage
    prepared = storage
    for key in Jutul.submodels_symbols(model)
        submodel = model[key]
        if submodel isa FacilityModel
            substorage = prepared[key]
            state = substorage[:state]
            state0 = substorage[:state0]
            T = eltype(state.FacilityCrossTermState.factor)
            state = merge(state, (FacilityCrossTermState =
                FacilityCrossTermState(submodel;
                    T = T, number_of_tracers = ntracers),))
            state0 = merge(state0, (FacilityCrossTermState =
                FacilityCrossTermState(submodel;
                    T = T, number_of_tracers = ntracers),))
            substorage = Jutul.JutulStorage(merge(Jutul.data(substorage),
                (; state, state0)))
            state_storage = Jutul.JutulStorage(merge(
                Jutul.data(prepared[:state]), NamedTuple{(key,)}((state,))))
            state0_storage = Jutul.JutulStorage(merge(
                Jutul.data(prepared[:state0]), NamedTuple{(key,)}((state0,))))
            prepared = Jutul.JutulStorage(merge(Jutul.data(prepared),
                NamedTuple{(key,)}((substorage,)),
                (; state = state_storage, state0 = state0_storage)))
        end
    end
    return prepared
end

function Adapt.adapt_structure(to, ct::ReservoirFromWellFlowCT)
    return ReservoirFromWellFlowCT(
        Adapt.adapt(to, ct.reservoir_cells),
        Adapt.adapt(to, ct.well_cells)
    )
end

function Adapt.adapt_structure(to, ct::ReservoirFromWellThermalCT)
    return ReservoirFromWellThermalCT(
        Adapt.adapt(to, ct.reservoir_cells),
        Adapt.adapt(to, ct.well_cells)
    )
end

function Adapt.adapt_structure(to, ct::WellFromFacilityFlowCT)
    return WellFromFacilityFlowCT(nothing, ct.facility_position)
end

function Adapt.adapt_structure(to, ct::WellFromFacilityThermalCT)
    return WellFromFacilityThermalCT(nothing, ct.facility_position)
end

function Adapt.adapt_structure(to, ct::FacilityFromWellTemperatureCT)
    return FacilityFromWellTemperatureCT(nothing, ct.facility_position)
end

function Adapt.adapt_structure(to, ct::FacilityFromWellEnthalpyCT)
    return FacilityFromWellEnthalpyCT(nothing, ct.facility_position)
end

function Adapt.adapt_structure(to, ct::FacilityFromWellBottomHolePressureCT)
    return FacilityFromWellBottomHolePressureCT(
        nothing, ct.facility_position)
end

function Adapt.adapt_structure(to, ct::FacilityFromSurfacePhaseRatesCT)
    return FacilityFromSurfacePhaseRatesCT(nothing, ct.facility_position)
end
