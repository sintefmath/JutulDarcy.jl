
struct SimpleCapillaryPressure{T, R} <: VectorVariables
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

function Jutul.line_plot_data(model::SimulationModel, cap::SimpleCapillaryPressure)
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

function Jutul.subvariable(p::SimpleCapillaryPressure, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(p.regions, c)
    return SimpleCapillaryPressure(p.pc, regions = regions)
end

degrees_of_freedom_per_entity(model, v::SimpleCapillaryPressure) = number_of_phases(model.system) - 1

@jutul_secondary function update_pc!(Δp, pc::SimpleCapillaryPressure, model, Saturations, ix)
    cap = pc.pc
    npc = size(Δp, 1)
    if npc == 1
        pcow = only(cap)
        @inbounds for c in ix
            reg = region(pc.regions, c)
            pcow_c = table_by_region(pcow, reg)
            sw = Saturations[1, c]
            Δp[1, c] = pcow_c(sw)
        end
    elseif npc == 2
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
                Δp[1, c] = -pcow_c(sw)
                Δp[2, c] = 0
            end
        else
            @inbounds for c in ix
                reg = region(pc.regions, c)
                pcow_c = table_by_region(pcow, reg)
                pcog_c = table_by_region(pcog, reg)
                sw = Saturations[1, c]
                sg = Saturations[3, c]
                Δp[1, c] = -pcow_c(sw)
                Δp[2, c] = pcog_c(sg)
            end
        end
    else
        error("Only implemented for two and three-phase flow.")
    end
end
