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

Adapt.adapt_structure(::Jutul.KernelAbstractionsContext, group::WellGroup) = group

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
