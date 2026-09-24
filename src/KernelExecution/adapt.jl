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
Adapt.@adapt_structure PVCDO
Adapt.@adapt_structure PVDG
Adapt.@adapt_structure PVTG
Adapt.@adapt_structure PVTGTable
Adapt.@adapt_structure PVTW
Adapt.@adapt_structure WATDENT

adapt_reservoir_scalar(ctx, value) = convert(
    Jutul.KernelExecution.ka_storage_eltype(ctx, typeof(value)), value)
Adapt.@adapt_structure FacilitySystem
Adapt.@adapt_structure KValueWrapper
Adapt.@adapt_structure PTViscosities
Adapt.@adapt_structure BrineCO2MixingDensities

Jutul.KernelExecution.ka_storage_eltype(ctx::Jutul.KernelAbstractionsContext,
    ::Type{BlackOilX{T}}) where T =
    BlackOilX{Jutul.KernelExecution.ka_storage_eltype(ctx, T)}

function Jutul.KernelExecution.ka_storage_eltype(
        ctx::Jutul.KernelAbstractionsContext,
        ::Type{MultiComponentFlash.FlashedMixture2Phase{T, A, E, R}}) where {T, A, E, R}
    F = Jutul.KernelExecution.ka_storage_eltype(ctx, T)
    V = Jutul.KernelExecution.ka_storage_eltype(ctx, A)
    K = Jutul.KernelExecution.ka_storage_eltype(ctx, E)
    return MultiComponentFlash.FlashedMixture2Phase{F, V, K, eltype(K)}
end

function Adapt.adapt_structure(ctx::Jutul.KernelAbstractionsContext,
        x::BlackOilX)
    T = Jutul.KernelExecution.ka_storage_eltype(ctx, typeof(x.val))
    return convert(BlackOilX{T}, x)
end

function Adapt.adapt_structure(ctx::Jutul.KernelAbstractionsContext,
        variable::BlackOilUnknown)
    F = Jutul.float_type(ctx)
    return BlackOilUnknown(
        dr_max = convert(F, variable.dr_max),
        ds_max = convert(F, variable.ds_max))
end

function Adapt.adapt_structure(to,
        variable::PressureTemperatureDependentVariable{T, R, N}) where {T, R, N}
    return PressureTemperatureDependentVariable(
        Adapt.adapt(to, variable.tab),
        Adapt.adapt(to, variable.regions),
        Val(N))
end

function adapt_compositional_eos(to, eos; float_type = missing)
    return MultiComponentFlash.make_eos_immutable(eos; float_type = float_type)
end

function adapt_compositional_eos(to, eos::MultiComponentFlash.KValuesEOS; float_type = missing)
    eos = MultiComponentFlash.make_eos_immutable(eos; float_type = float_type)
    evaluator = Adapt.adapt(to, eos.K_values_evaluator)
    return MultiComponentFlash.KValuesEOS(evaluator, eos.mixture;
        volume_shift = eos.volume_shift)
end

function Adapt.adapt_structure(to,
        system::MultiPhaseCompositionalSystemLV{E, T, O, R, N, C, Ref}) where {
        E, T, O, R, N, C, Ref}
    if to isa Jutul.JutulContext
        float_type = Jutul.float_type(to)
    else
        float_type = missing
    end
    eos = adapt_compositional_eos(to, system.equation_of_state, float_type = float_type)
    phases = Adapt.adapt(to, system.phases)
    rho_ref = Adapt.adapt(to, system.rho_ref)
    return MultiPhaseCompositionalSystemLV{
        typeof(eos), typeof(phases), O, typeof(rho_ref), N, Nothing, Ref}(
        phases, nothing, eos, rho_ref, system.reference_phase_index)
end

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

Adapt.adapt_structure(::Jutul.KernelAbstractionsContext, group::WellGroup) = group

function Adapt.adapt_structure(to, well::SimpleWell)
    return SimpleWell(
        Adapt.adapt(to, well.perforations),
        Adapt.adapt(to, well.surface),
        nothing,
        well.explicit_dp,
        well.multiwell
    )
end

function Adapt.adapt_structure(to,
        system::StandardBlackOilSystem{D, V, W, R, F, T, P, Num, Ref}) where {D, V, W, R, F, T, P, Num, Ref}
    rs_max = Adapt.adapt(to, system.rs_max)
    rv_max = Adapt.adapt(to, system.rv_max)
    rho_ref = Adapt.adapt(to, system.rho_ref)
    phase_indices = Adapt.adapt(to, system.phase_indices)
    phases = Adapt.adapt(to, system.phases)
    if to isa Jutul.KernelAbstractionsContext
        Num_t = Jutul.float_type(to)
    else
        Num_t = Num
    end
    rho_ref = map(x -> convert(Num_t, x), rho_ref)
    return StandardBlackOilSystem{
        typeof(rs_max), typeof(rv_max), W, typeof(rho_ref), F,
        typeof(phase_indices), typeof(phases), Num_t, Ref
    }(rs_max, rv_max, rho_ref, phase_indices, phases,
        system.saturated_chop, system.keep_bubble_flag,
        convert(Num_t, system.rs_eps), convert(Num_t, system.rv_eps),
        convert(Num_t, system.s_eps),
        system.reference_phase_index)
end

function Adapt.adapt_structure(ctx::Jutul.KernelAbstractionsContext,
        system::ImmiscibleSystem{T, F, Ref}) where {T, F, Ref}
    rho_ref = map(x -> adapt_reservoir_scalar(ctx, x), system.rho_ref)
    return ImmiscibleSystem{T, typeof(rho_ref), Ref}(
        system.phases, rho_ref, system.reference_phase_index)
end

function Adapt.adapt_structure(to,
        variable::ReservoirRelativePermeabilities{Scaling, ph}) where {Scaling, ph}
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
    return ReservoirRelativePermeabilities{
        typeof(scaling), ph, typeof(krw), typeof(krow), typeof(krog), typeof(krg),
        typeof(regions), typeof(hysteresis_w), typeof(hysteresis_ow),
        typeof(hysteresis_og), typeof(hysteresis_g), typeof(method)
    }(krw, krow, krog, krg, regions, hysteresis_w, hysteresis_ow,
        hysteresis_og, hysteresis_g, scaling,
        variable.hysteresis_s_threshold, variable.hysteresis_s_eps, method)
end

function Adapt.adapt_structure(to, variable::PhaseRelativePermeability)
    return PhaseRelativePermeability(
        Adapt.adapt(to, variable.k), variable.label,
        variable.connate, variable.critical, variable.s_max,
        variable.k_max, variable.input_s_max)
end

function Adapt.adapt_structure(ctx::Jutul.KernelAbstractionsContext,
        variable::PhaseRelativePermeability)
    return PhaseRelativePermeability(
        Adapt.adapt(ctx, variable.k), variable.label,
        adapt_reservoir_scalar(ctx, variable.connate),
        adapt_reservoir_scalar(ctx, variable.critical),
        adapt_reservoir_scalar(ctx, variable.s_max),
        adapt_reservoir_scalar(ctx, variable.k_max),
        adapt_reservoir_scalar(ctx, variable.input_s_max))
end

function Adapt.adapt_structure(ctx::Jutul.KernelAbstractionsContext,
        table::ConstMuBTable)
    return ConstMuBTable(
        adapt_reservoir_scalar(ctx, table.p_ref),
        adapt_reservoir_scalar(ctx, table.b_ref),
        adapt_reservoir_scalar(ctx, table.b_c),
        adapt_reservoir_scalar(ctx, table.mu_ref),
        adapt_reservoir_scalar(ctx, table.mu_c))
end

function Adapt.adapt_structure(ctx::Jutul.KernelAbstractionsContext,
        table::WATDENT{N}) where N
    tab = map(table.tab) do record
        (T = adapt_reservoir_scalar(ctx, record.T),
            c1 = adapt_reservoir_scalar(ctx, record.c1),
            c2 = adapt_reservoir_scalar(ctx, record.c2))
    end
    return WATDENT{N, typeof(first(tab))}(tab)
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

function Adapt.adapt_structure(ctx::Jutul.KernelAbstractionsContext,
        variable::LinearlyCompressiblePoreVolume)
    reference_pressure = map(x -> adapt_reservoir_scalar(ctx, x),
        variable.reference_pressure)
    expansion = map(x -> adapt_reservoir_scalar(ctx, x),
        variable.expansion)
    return LinearlyCompressiblePoreVolume(
        Adapt.adapt(ctx, reference_pressure),
        Adapt.adapt(ctx, expansion),
        Adapt.adapt(ctx, variable.regions),
        Val(:assembled))
end

function Adapt.adapt_structure(to, mask::PerforationMask)
    values = Adapt.adapt(to, mask.values)
    return PerforationMask(values, Val(:adapted))
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
