abstract type AbstractCapillaryPressure <: VectorVariables end

degrees_of_freedom_per_entity(model, v::AbstractCapillaryPressure) = number_of_phases(model.system) - 1

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
        is_tup_tup = first(pc) isa Tuple
        if isnothing(regions)
            @assert !is_tup_tup || all(x -> length(x) == 1, pc)
        end
        pc = map(x -> region_wrap(x, regions), pc)
        pc = tuple(pc...)
        return new{typeof(pc), T}(pc, regions)
    end
end

function Jutul.subvariable(p::SimpleCapillaryPressure, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(p.regions, c)
    return SimpleCapillaryPressure(p.pc, regions = regions)
end


@jutul_secondary function update_pc!(Δp, pc::SimpleCapillaryPressure, model, Saturations, ix)
    cap = pc.pc
    npc = size(Δp, 1)
    reference_ph = get_reference_phase_index(model.system)
    if npc == 1
        pcow = only(cap)
        @assert reference_ph == 1
        @inbounds for c in ix
            reg = region(pc.regions, c)
            pcow_c = table_by_region(pcow, reg)
            sw = Saturations[1, c]
            Δp[1, c] = pcow_c(sw)
        end
    elseif npc == 2
        @assert reference_ph == 2
        pcow, pcog = cap
        if isnothing(pcow)
            @inbounds for c in ix
                reg = region(pc.regions, c)
                pcog_c = table_by_region(pcog, reg)
                sg = Saturations[3, c]
                Δp[1, c] = 0
                Δp[2, c] = pcog_c(sg)
            end
        elseif isnothing(pcog)
            @inbounds for c in ix
                reg = region(pc.regions, c)
                pcow_c = table_by_region(pcow, reg)
                sw = Saturations[1, c]
                Δp[1, c] = pcow_c(sw)
                Δp[2, c] = 0
            end
        else
            @inbounds for c in ix
                reg = region(pc.regions, c)
                pcow_c = table_by_region(pcow, reg)
                pcog_c = table_by_region(pcog, reg)
                sw = Saturations[1, c]
                sg = Saturations[3, c]
                # Note: Negative sign already taken care of in input
                Δp[1, c] = pcow_c(sw)
                Δp[2, c] = pcog_c(sg)
            end
        end
    else
        error("Only implemented for two and three-phase flow.")
    end
end

struct ScaledCapillaryPressure{T, R, S} <: AbstractCapillaryPressure
    pc::T
    regions::R
    scaling::S
    function ScaledCapillaryPressure(pc::C, scaling::S; regions::T = nothing) where {C, S, T}
        is_tup_tup = first(pc) isa Tuple
        if isnothing(regions)
            @assert !is_tup_tup || all(x -> length(x) == 1, pc)
        end
        pc = map(x -> region_wrap(x, regions), pc)
        pc = tuple(pc...)
        @assert length(pc) == size(scaling, 1)
        return new{typeof(pc), T, S}(pc, regions, scaling)
    end
end

function Jutul.subvariable(p::ScaledCapillaryPressure, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(p.regions, c)
    scaling = Jutul.partition_variable_slice(p.scaling, c)
    return ScaledCapillaryPressure(p.pc, scaling, regions = regions)
end

@jutul_secondary function update_pc!(Δp, pc::ScaledCapillaryPressure, model, Saturations, ix)
    cap = pc.pc
    npc = size(Δp, 1)
    scale = pc.scaling
    if npc == 1
        pcow = only(cap)
        @inbounds for c in ix
            reg = region(pc.regions, c)
            pcow_c = table_by_region(pcow, reg)
            sw = Saturations[1, c]
            Δp[1, c] = scale[1, c]*pcow_c(sw)
        end
    elseif npc == 2
        pcow, pcog = cap
        if isnothing(pcow)
            @inbounds for c in ix
                reg = region(pc.regions, c)
                pcog_c = table_by_region(pcog, reg)
                sg = Saturations[3, c]
                Δp[1, c] = 0
                Δp[2, c] = scale[2, c]*pcog_c(sg)
            end
        elseif isnothing(pcog)
            @inbounds for c in ix
                reg = region(pc.regions, c)
                pcow_c = table_by_region(pcow, reg)
                sw = Saturations[1, c]
                Δp[1, c] = -scale[1, c]*pcow_c(sw)
                Δp[2, c] = 0
            end
        else
            @inbounds for c in ix
                reg = region(pc.regions, c)
                pcow_c = table_by_region(pcow, reg)
                pcog_c = table_by_region(pcog, reg)
                sw = Saturations[1, c]
                sg = Saturations[3, c]
                Δp[1, c] = -scale[1, c]*pcow_c(sw)
                Δp[2, c] = scale[2, c]*pcog_c(sg)
            end
        end
    else
        error("Only implemented for two and three-phase flow.")
    end
end
