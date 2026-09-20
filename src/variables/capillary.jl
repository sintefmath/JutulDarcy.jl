abstract type AbstractCapillaryPressure <: VectorVariables end

degrees_of_freedom_per_entity(model, v::AbstractCapillaryPressure) = number_of_phases(model.system) - 1

# Use the same interpolation representation as a supplied table, including its
# lookup mode and storage type, for an omitted capillary pressure curve.
function zero_capillary_interpolator(table::Jutul.LinearInterpolant)
    return Jutul.LinearInterpolant(table.X, zero.(table.F), table.lookup)
end

zero_capillary_interpolator(::Nothing) =
    get_1d_interpolator([0.0, 1.0], [0.0, 0.0])

function capillary_template(pc)
    for pair in pc
        tables = pair isa Union{Tuple, AbstractVector} ? pair : (pair,)
        for table in tables
            if table isa Jutul.LinearInterpolant
                return table
            end
        end
    end
    return nothing
end

function normalize_capillary_pair(pair, fallback, nregions)
    if pair isa Union{Tuple, AbstractVector}
        any(isnothing, pair) || return pair
        template = capillary_template((pair,))
        zero_table = zero_capillary_interpolator(isnothing(template) ? fallback : template)
        return map(table -> isnothing(table) ? zero_table : table, pair)
    else
        if isnothing(pair)
            zero_table = zero_capillary_interpolator(fallback)
            if nregions == 1
                return zero_table
            elseif nregions < 20
                return ntuple(_ -> zero_table, nregions)
            else
                return fill(zero_table, nregions)
            end
        end
        return pair
    end
end

function prepare_capillary_tables(pc, regions)
    template = capillary_template(pc)
    nregions = regions isa AbstractArray ? maximum(regions) : 1
    return tuple((region_wrap(normalize_capillary_pair(pair, template, nregions), regions) for pair in pc)...)
end

function Jutul.line_plot_data(model::SimulationModel, cap::AbstractCapillaryPressure)
    npc = number_of_phases(model.system)-1
    phases = phase_names(model.system)
    nreg = length(cap.pc[1])
    data = Matrix{Any}(undef, 1, nreg)
    for reg in 1:nreg
        x = []
        y = []
        labels = []
        for i in 1:npc
            pc = cap.pc[i]
            (; X, F) = pc[reg]
            push!(x, X[2:end-1])
            push!(y, F[2:end-1]./1e5)
            prev = phases[i]
            next = phases[i+1]
            push!(labels, "$prev-$next")
        end
        data[reg] = JutulLinePlotData(x, y, title = "Capillary pressure", xlabel = "Saturation", ylabel = "Pc [bar]", labels = labels)
    end
    return data
end

struct SimpleCapillaryPressure{T, R} <: AbstractCapillaryPressure
    pc::T
    regions::R
    function SimpleCapillaryPressure(pc::C; regions::T = nothing) where {C, T}
        if isnothing(regions)
            @assert all(x -> !(x isa Union{Tuple, AbstractVector}) || length(x) == 1, pc)
        end
        pc = prepare_capillary_tables(pc, regions)
        return new{typeof(pc), T}(pc, regions)
    end
    function SimpleCapillaryPressure(pc::T, regions::R,
            ::Val{:assembled}) where {T, R}
        return new{T, R}(pc, regions)
    end
end

function SimpleCapillaryPressure(pc::Jutul.LinearInterpolant; kwarg...)
    return SimpleCapillaryPressure((pc, ); kwarg...)
end

function Jutul.subvariable(p::SimpleCapillaryPressure, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(p.regions, c)
    return SimpleCapillaryPressure(p.pc, regions = regions)
end


@jutul_secondary function update_pc!(Δp, pc::SimpleCapillaryPressure, model, Saturations, ix)
    return update_capillary_pressure!(Δp, pc, model, Saturations, nothing, ix)
end

struct CapillaryPressureScaling <: Jutul.VectorVariables end

Jutul.degrees_of_freedom_per_entity(model, v::CapillaryPressureScaling) = number_of_phases(model.system) - 1
Jutul.default_value(::CapillaryPressureScaling) = 1.0

function Jutul.default_parameter_values(data_domain, model, param::CapillaryPressureScaling, symb)
    N = Jutul.degrees_of_freedom_per_entity(model, param)
    M = number_of_cells(model.domain)
    if haskey(data_domain, :capillary_pressure_scaling)
        out = data_domain[:capillary_pressure_scaling, Cells()]
    else
        v = Jutul.default_value(param)
        out = fill(v, N, M)
    end
    size(out) == (N, M) || error("Capillary pressure scaling must have size (N, number of cells), where N is the number of non-reference phases.")
    return out
end

struct ScaledCapillaryPressure{T, R} <: AbstractCapillaryPressure
    pc::T
    regions::R
    function ScaledCapillaryPressure(pc::C; regions::T = nothing) where {C, T}
        if isnothing(regions)
            @assert all(x -> !(x isa Union{Tuple, AbstractVector}) || length(x) == 1, pc)
        end
        pc = prepare_capillary_tables(pc, regions)
        return new{typeof(pc), T}(pc, regions)
    end
    function ScaledCapillaryPressure(pc::T, regions::R,
            ::Val{:assembled}) where {T, R}
        return new{T, R}(pc, regions)
    end
end

function Jutul.subvariable(p::ScaledCapillaryPressure, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(p.regions, c)
    return ScaledCapillaryPressure(p.pc, regions = regions)
end

@jutul_secondary function update_pc!(Δp, pc::ScaledCapillaryPressure, model, Saturations, CapillaryPressureScaling, ix)
    return update_capillary_pressure!(Δp, pc, model, Saturations, CapillaryPressureScaling, ix)
end

@inline capillary_scale(::Nothing, phase, cell) = 1
@inline capillary_scale(scale, phase, cell) = scale[phase, cell]

function update_capillary_pressure!(Δp, pc, model, Saturations, scale, ix)
    cap = pc.pc
    npc = size(Δp, 1)
    reference_ph = get_reference_phase_index(model.system)
    if npc == 1
        if reference_ph == 1
            w = 2
        else
            w = 1
        end
        pcow = cap[1]
        @inbounds for c in ix
            reg = region(pc.regions, c)
            sw = Saturations[w, c]
            Δp[1, c] = capillary_scale(scale, 1, c)*evaluate_table_by_region(pcow, reg, sw)
        end
    elseif npc == 2
        if reference_ph == 1
            w, g = 2, 3
        elseif reference_ph == 2
            w, g = 1, 3
        else
            @assert reference_ph == 3
            w, g = 1, 2
        end
        pcow, pcog = cap
        @inbounds for c in ix
            reg = region(pc.regions, c)
            sw = Saturations[w, c]
            sg = Saturations[g, c]
            # The input interpolators already include the capillary pressure sign.
            Δp[1, c] = capillary_scale(scale, 1, c)*evaluate_table_by_region(pcow, reg, sw)
            Δp[2, c] = capillary_scale(scale, 2, c)*evaluate_table_by_region(pcog, reg, sg)
        end
    else
        error("Only implemented for two and three-phase flow.")
    end
end


"""
    pc = brooks_corey_pc(s;
        p_entry = 2e5,
        n = 2.0,
        p_max = Inf,
        p_min = -Inf,
        residual = 0.0,
        residual_total = residual
    )

Compute the capillary pressure for a given saturation using the Brooks-Corey model.

# Arguments
- `s::Real`: Saturation

# Keyword Arguments
- `p_entry::Real`: Entry pressure
- `n::Real`: Corey exponent
- `p_max::Real`: Maximum capillary pressure
- `residual::Real`: Residual saturation
- `residual_total::Real`: Total residual saturation
"""
function brooks_corey_pc(s;
        p_entry = 2e5,
        n = 2.0,
        p_max = Inf,
        residual = 0.0,
        residual_total = residual
    )
    @assert s <= 1.0
    @assert s >= 0.0
    @assert residual <= 1.0
    @assert residual >= 0.0
    @assert residual_total <= 1.0
    @assert residual_total >= 0.0

    @assert residual <= residual_total
    @assert isfinite(p_entry)
    return brooks_corey_pc(s, p_entry, n, residual, residual_total, p_max)
end

function brooks_corey_pc(s, p_e, n, residual, residual_total, p_max = Inf)
    s_norm = normalized_saturation(s, residual, residual_total)
    return min(p_e*s_norm^(-1.0/n), p_max)
end
