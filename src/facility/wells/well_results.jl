function Base.show(io::IO, t::MIME"text/plain", wr::WellResults)
    n = length(wr.time)
    wells = keys(wr.wells)
    if length(wells) == 0
        outputs = Symbol[]
    else
        outputs = keys(wr.wells[first(wells)])
    end
    println(io, "WellResults with $n entries:")
    println(io, "  Wells: $(join(wells, ", ",))")
    println(io, "  Outputs: $(join(outputs, ", ",))")

    println(io, "\n  Properties:")
    println(io, "    .start_date = $(wr.start_date)")
    println(io, "    .time = $(typeof(wr.time))")
    println(io, "    .wells = $(typeof(wr.wells))")
end

function Base.getindex(wr::WellResults, well::Symbol)
    # Get all results for well
    return wr.wells[well]
end

function Base.getindex(wr::WellResults, well::Symbol, output::Symbol)
    # Single well and single output
    return wr.wells[well][output]
end

function Base.getindex(wr::WellResults, well::Symbol, outputs::Base.AbstractVecOrTuple)
    # Single well and multiple outputs
    n = length(wr.time)
    m = length(outputs)
    ret = zeros(n, m)
    for (i, output) in enumerate(outputs)
        ret[:, i] .= wr[well, output]
    end
    return ret
end

function Base.getindex(wr::WellResults, wells::Base.AbstractVecOrTuple, output::Symbol)
    # Multiple wells, single output
    n = length(wr.time)
    m = length(wells)
    ret = zeros(n, m)
    for (i, well) in enumerate(wells)
        ret[:, i] .= wr[well, output]
    end
    return ret
end

function Base.getindex(wr::WellResults, wells::Base.AbstractVecOrTuple, outputs::Base.AbstractVecOrTuple)
    # Multiple wells, multiple outputs
    return map(w -> wr[w, outputs], wells)
end

##



function well_result_dates(wr::WellResults)
    n = length(wr.time)
    if isnothing(wr.start_date)
        header = (["time"], ["days"])
        data = zeros(n, 1)
        data .= wr.time./si_unit(:day)
    else
        header = (["time", "date"], ["days", "-"])
        dates = @. DateTime(wr.start_date) + Second(wr.time)
        data = hcat(wr.time, dates)
    end
    return (header = header, data = data, title = "Time")
end

function well_result_output_legend_table(wr::WellResults, outputs = missing)
    if ismissing(outputs)
        outputs = collect(keys(wr.wells[first(keys(wr.wells))]))
    end
    sort!(outputs, by = x -> "$x")
    data = Matrix{Any}(undef, length(outputs), 4)
    for (i, output) in enumerate(outputs)
        info = JutulDarcy.well_target_information(output)
        if ismissing(info)
            data[i, 1] = output
            data[i, 2] = "?"
            data[i, 3] = "?"
            data[i, 4] = "?"
        else
            data[i, 1] = info.symbol
            data[i, 2] = info.description
            # data[i, 3] = info.explanation
            data[i, 3] = info.unit_label
            if info.is_rate
                ut = "$(info.unit_type) per time"
            else
                ut = "$(info.unit_type)"
            end
            data[i, 4] = ut
        end
    end

    header = ["Label", "Description", "Unit", "Type of quantity"]
    return (header = header, data = data, title = "Legend")
end

function well_result_table(wr::WellResults, well::Symbol, outputs = collect(keys(wr.wells[well])))
    time_tbl = well_result_dates(wr)
    r = map(x -> wr[well, x], outputs)
    data = hcat(time_tbl.data, r...)
    header = time_tbl.header
    for output in outputs
        descr = JutulDarcy.well_target_information(output)
        if ismissing(descr)
            u = "-"
        else
            u = descr.unit_label
        end
        push!(header[1], "$output")
        push!(header[2], u)
    end

    return (header = header, data = data, title = "$well result")
end

function well_result_table(wr::WellResults, well::Symbol, output::Symbol)
    well_result_table(wr, [well], output)
end

function well_result_table(wr::WellResults, wells::Base.AbstractVecOrTuple, output::Symbol)
    time_tbl = well_result_dates(wr)
    r = map(well -> wr[well, output], wells)
    data = hcat(time_tbl.data, r...)
    header = time_tbl.header
    descr = JutulDarcy.well_target_information(output)
    if ismissing(descr)
        t = "$output"
        u = "-"
    else
        t = "$(descr.description) ($output)"
        u = descr.unit_label
    end
    for well in wells
        push!(header[1], "$well")
        push!(header[2], "$u")
    end
    return (header = header, data = data, title = t)
end

"""
    print_well_result_table(wr::WellResults, wells)
    print_well_result_table(wr::WellResults, wells, outputs)

Print summary tables that show the well responses.
"""
function print_well_result_table(wr::WellResults, wells, arg...; kwarg...)
    # Wrapper to call with stdout
    print_well_result_table(stdout, wr, wells, arg...; kwarg...)
end

function print_well_result_table(io::IO, wr::WellResults, wells, arg...; legend = true, kwarg...)
    # Main caller
    function print_well(w)
        tbl = well_result_table(wr, w, arg...; kwarg...)
        Jutul.PrettyTables.pretty_table(io, tbl.data, header = tbl.header, title = tbl.title)
    end
    wells_is_iterable = wells isa Base.AbstractVecOrTuple
    outputs_is_iterable = length(arg) == 0 || first(arg) isa Base.AbstractVecOrTuple

    if outputs_is_iterable && legend
        tab = well_result_output_legend_table(wr, arg...)
        Jutul.PrettyTables.pretty_table(io, tab.data,
            header = tab.header,
            title = tab.title,
            alignment = :l,
        )
    end
    if wells_is_iterable && outputs_is_iterable
        # Both args are iterable, loop over wells
        for well in wells
            print_well(well)
        end
    else
        print_well(wells)
    end
end

function (wr::WellResults)(wells = collect(keys(wr.wells)), outputs = missing; kwarg...)
    if wells isa Symbol
        wells = [wells]
    elseif wells isa Colon
        wells = collect(keys(wr.wells))
    end
    if length(wells) > 0
        if ismissing(outputs)
            outputs = collect(keys(wr.wells[first(wells)]))
            sort!(outputs, by = x -> "$x")
        elseif outputs isa Symbol
            outputs = [outputs]
        end
        print_well_result_table(wr, wells, outputs; kwarg...)
    end
end
