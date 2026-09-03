
function JutulDarcy.plot_summary_impl(v::Union{AbstractVector, Tuple}; kwarg...)
    return JutulDarcy.plot_summary_impl(v...; kwarg...)
end

function JutulDarcy.plot_summary_impl(arg...;
        names = missing,
        unit_system = "Metric",
        linewidth = 2.0,
        markersize = 0.0,
        plots = ["FIELD:FPR"],
        extra_field = String[],
        extra_well = String[],
        events = missing,
        linecolor = :black,
        selectors = true,
        plot_type = :lines,
        cols = 1,
        alpha = 1.0,
        nxticks = 15,
        legend = true,
        hideaxes = false,
        rows = ceil(length(plots)/cols) |> Int,
        colormap = missing,
        kwarg...
    )
    if length(arg) == 1 && only(arg) isa AbstractVector
        arg = only(arg)
    end
    if ismissing(names)
        names = ["Summary $i" for i in 1:length(arg)]
    end
    function split_name(inp::String)
        sep = ':'
        if in(sep, inp)
            well_or_fld, name = split(inp, ':')
        else
            well_or_fld = "FIELD"
            name = inp
        end
        return (well_or_fld, name)
    end
    extra_well_internal = String[]
    extra_field_internal = String[]
    for k in plots
        well_or_fld, name = split_name(k)
        dest = well_or_fld == "FIELD" ? extra_field_internal : extra_well_internal
        if !(name in dest)
            push!(dest, name)
        end
    end
    for (prefix, dest) in zip(["F", "W"], [extra_field_internal, extra_well_internal])
        push!(dest, "$(prefix)WPR,$(prefix)OPR,$(prefix)GPR")
        push!(dest, "$(prefix)WPT,$(prefix)OPT,$(prefix)GPT")
        push!(dest, "$(prefix)WIT,$(prefix)WPT")
        push!(dest, "$(prefix)OIP,$(prefix)OPT")
        push!(dest, "$(prefix)WIP,$(prefix)WPT")
    end

    get_summary(r::JutulDarcy.ReservoirSimResult) = r.summary
    get_summary(x) = x
    summaries = [get_summary(s) for s in arg]
    nsmry = length(names)
    nsmry == length(summaries) || error("Number of names ($nsmry) must match number of summaries ($(length(summaries))).")
    for (i, s) in enumerate(summaries)
        summaries[i] = copy(s)
        if !haskey(s, "UNIT_SYSTEM")
            println("Warning: Summary $i has no UNIT_SYSTEM key, assuming metric.")
            s["UNIT_SYSTEM"] = "metric"
        end
        if !haskey(s["VALUES"], "FIELD")
            println("Warning: Summary $i has no FIELD values, adding empty.")
            s["VALUES"]["FIELD"] = Dict{String, Any}()
        end
    end
    if ismissing(colormap)
        wong_c = Makie.wong_colors()
        if nsmry <= length(wong_c)
            colormap = to_colormap(wong_c)
        else
            colormap = to_colormap(:tab20)
        end
    elseif isa(colormap, Symbol)
        colormap = to_colormap(colormap)
    end
    summary_sample = summaries[1]
    # The source of the data
    well_names = collect(keys(summary_sample["VALUES"]["WELLS"]))

    source_keys = copy(well_names)
    for i in 2:nsmry
        # Find intersection of well names
        smry = summaries[i]
        wnames_i = collect(keys(smry["VALUES"]["WELLS"]))
        source_keys = intersect(source_keys, wnames_i)
    end

    sort!(source_keys)
    pushfirst!(source_keys, "FIELD")
    # Field types
    field_quantity_keys = collect(keys(summary_sample["VALUES"]["FIELD"]))
    sort!(field_quantity_keys)
    pushfirst!(field_quantity_keys, "NONE")

    lookup = JutulDarcy.summary_key_lookup()
    function get_well_quantity_keys(wname)
        wk = collect(keys(summary_sample["VALUES"]["WELLS"][wname]))
        for k in cat(extra_well, extra_well_internal; dims = 1)
            if !(k in wk)
                push!(wk, k)
            end
        end
        return wk
    end

    function get_quantity_options(kind)
        if kind == "FIELD"
            opts = copy(field_quantity_keys)
            for k in cat(extra_field, extra_field_internal; dims = 1)
                if !(k in opts)
                    push!(opts, k)
                end
            end
        else
            opts = get_well_quantity_keys(kind)
        end
        return opts
    end

    function time_data(idx; maybe_datetime = false)
        smry = summaries[idx]
        t = smry["TIME"].seconds
        if maybe_datetime && !isnothing(start_date)
            out = start_dates[idx] .+ @. Microsecond(ceil(t*1e6))
        else
            out = t./si_unit(:day)
        end
        return out
    end

    function plot_data(kind, valtype, idx, info, units)
        smry = summaries[idx]
        t = smry["TIME"].seconds
        if valtype == "NONE"
            v = zeros(length(t))
        elseif kind == "FIELD"
            v = get(smry["VALUES"]["FIELD"], valtype, missing)
        else
            v = get(smry["VALUES"]["WELLS"][kind], valtype, missing)
        end
        if ismissing(v)
            println("Warning: No data for $kind:$valtype in $(names[idx])")
            return fill(NaN, length(t))
        end
        v = copy(v)
        from_sys = JutulDarcy.GeoEnergyIO.InputParser.DeckUnitSystem(Symbol(lowercase(smry["UNIT_SYSTEM"])))
        to_sys = JutulDarcy.GeoEnergyIO.InputParser.DeckUnitSystem(Symbol(lowercase(units)))
        systems = (to = to_sys, from = from_sys)
        if !ismissing(info) && v isa AbstractArray
            JutulDarcy.GeoEnergyIO.InputParser.swap_unit_system!(v, systems, info.unit_type)
        end
        return v
    end

    # Set up time labels
    start_dates = map(s -> s["TIME"].start_date, summaries)
    filtered_start_dates = copy(start_dates)
    unique!(filter!(x -> !isnothing(x) && !ismissing(x), filtered_start_dates))
    has_missing = length(filtered_start_dates) < length(start_dates)
    if length(filtered_start_dates) == 0
        start_date = nothing
    else
        start_date = first(filtered_start_dates)
        if has_missing
            if length(filtered_start_dates) > 1
                println("Note: Multiple distinct start dates ($start_dates) found in summaries, using first entry as the base for values without dates.")
            end
            for (i, s) in enumerate(summaries)
                if isnothing(s["TIME"].start_date) || ismissing(s["TIME"].start_date)
                    start_dates[i] = start_date
                end
            end
        end
    end
    # Max time - in days
    plots = map(String, plots)
    fig = Figure(size = (1200, 800))
    top_menu_grid = GridLayout(fig[2, 1])
    row_menu, = label_menu(top_menu_grid[1, 1], collect(range(1, max(rows, 5))), "Number of rows", default = string(rows))
    col_menu, = label_menu(top_menu_grid[1, 2], collect(range(1, max(cols, 5))), "Number of columns", default = string(cols))

    lineoptions = collect(range(0.0, 8.0, 17))
    push!(lineoptions, linewidth)
    unique!(lineoptions)
    sort!(lineoptions)

    markeroptions = [0.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
    push!(markeroptions, markersize)
    unique!(markeroptions)
    sort!(markeroptions)

    linewidth_menu, = label_menu(top_menu_grid[1, 3], lineoptions, "Line width", default = string(linewidth))
    menu_idx = 4
    if plot_type == :scatter || plot_type == :scatterlines
        markersize_menu, = label_menu(top_menu_grid[1, menu_idx], markeroptions, "Marker size", default = string(markersize))
        ms_sel = markersize_menu.selection
        menu_idx += 1
    else
        ms_sel = missing
    end

    unit_menu, = label_menu(top_menu_grid[1, menu_idx], ["Metric", "SI", "Field"], "Unit system"; default = unit_system)

    plot_layout = nothing
    plot_boxes = []
    legends = Dict{Int, Legend}()

    function plot_string(source::String, type::String)
        return "$source:$type"
    end

    function plot_string(source::String, ::Nothing)
        return plot_string(source, "NONE")
    end

    function float_fmt(x, u)
        if x > 1e5
            numstr = @sprintf("%g", x)
        else
            numstr = "$(round(x, sigdigits = 4))"
        end
        return "$numstr\n$u"
    end

    function update_plots(idx = eachindex(plots))
        for i in idx
            ax = plot_boxes[i].ax
            empty!(ax)
            well_or_fld, name = split_name(plots[i])

            plot_names = split(name, ',')
            ystr = ""
            tstr = ""
            nplts = length(plot_names)
            nlines = nplts*nsmry
            for (pn_idx, pname) in enumerate(plot_names)
                info = get(lookup, pname, missing)
                units = unit_menu.selection[]

                for (smry_no, smry_name) in enumerate(names)
                    t = time_data(smry_no, maybe_datetime = true)
                    v = plot_data(well_or_fld, pname, smry_no, info, units)
                    if !ismissing(info) && name != "NONE" && v isa AbstractArray
                        # ax.title[] = "$(name) $(info.legend)"
                        _, u = well_unit_conversion(units, "", info)
                        if info.is_rate
                            if lowercase(string(summaries[smry_no]["UNIT_SYSTEM"])) == "si"
                                v = v.*si_unit(:day)
                            end
                            u *= "/day"
                        end
                        if smry_no == 1
                            # Only add if we are the first summary
                            ttxt = "$(pname) $(info.legend)"
                            if pn_idx == 1
                                ystr = "$u"
                                tstr = ttxt
                            else
                                ystr = "$ystr, $u"
                                tstr = "$tstr, $ttxt"
                            end
                        end
                    end
                    lw_sel = linewidth_menu.selection
                    if nplts == 1
                        if nsmry == 1
                            arg = (color = linecolor, )
                        else
                            maxsmry = max(nsmry, 2)
                            arg = (color = smry_no, colorrange = (1, maxsmry), colormap = colormap[1:min(length(colormap), maxsmry)])
                        end
                    else
                        # Multiple summaries and multiple plots
                        # Use plot type to color
                        # Line type from summary
                        linestyles = (:solid, :dash, :dot, :dashdot, :dashdotdot)
                        ls = linestyles[mod1(smry_no, length(linestyles))]
                        maxplt = max(nplts, 2)
                        arg = (color = pn_idx, colorrange = (1, maxplt), colormap = colormap[1:min(length(colormap), maxplt)], linestyle = ls)
                    end
                    if nplts == 1
                        lbl = smry_name
                    elseif nsmry == 1
                        lbl = pname
                    else
                        lbl = "$(smry_name):$(pname)"
                    end
                    if all(x -> !isfinite(x) || x == 0.0, v)
                        lbl = lbl * " (no data)"
                    end
                    plot_arg = (
                        label = lbl,
                        alpha = alpha,
                        transparency = alpha < 1.0,
                        inspectable = false,
                        kwarg...
                    )
                    if plot_type == :lines
                        lines!(ax, t, v;
                            linewidth = lw_sel,
                            plot_arg...,
                            arg...
                        )
                    elseif plot_type == :stairs
                        stairs!(ax, t, v;
                            linewidth = lw_sel,
                            plot_arg...,
                            arg...
                        )
                    elseif plot_type == :scatter
                        scatter!(ax, t, v;
                            markersize = ms_sel,
                            plot_arg...,
                            arg...
                        )
                    elseif plot_type == :scatterlines
                        scatterlines!(ax, t, v;
                            linewidth = lw_sel,
                            markersize = ms_sel,
                            plot_arg...,
                            arg...
                        )
                    else
                        error("Unknown plot type: $plot_type")
                    end
                    if !ismissing(events)
                        plot_events!(ax, events, well_or_fld, start_date, t, v)
                    end
                end
                # ax.ylabel[] = ystr
                ax.ytickformat[] = values -> [float_fmt(value, ystr) for value in values]
                ax.title[] = tstr
            end
            if nlines > 1 && legend
                l = axislegend(ax, position = :lt)
                if haskey(legends, i)
                    delete!(legends[i])
                end
                legends[i] = l
            end
        end
    end


    function update_field_or_well_selection(m1, m2, new_selection, idx)
        m2.options[] = get_quantity_options(new_selection)
        plots[idx] = plot_string(new_selection, m2.selection[])
        update_plots(idx)
    end

    function update_response_selection(m1, m2, new_selection, idx)
        plots[idx] = plot_string(m1.selection[], new_selection)
        update_plots(idx)
    end

    function update_menu_layout()
        nrows = row_menu.selection[]
        ncols = col_menu.selection[]
        plot_layout = GridLayout(nrows, ncols)
        fig.layout[1, 1] = plot_layout
        len = nrows*ncols
        nplots = length(plots)
        for _ in 1:(len - nplots)
            push!(plots, "NONE")
        end
        resize!(plots, len)

        for el in plot_boxes
            delete!(el.menu1)
            delete!(el.menu2)
            delete!(el.label1)
            delete!(el.label2)
            delete!(el.ax)
        end
        for el in values(legends)
            delete!(el)
        end
        plot_boxes = Matrix{Any}(undef, nrows, ncols)
        plot_idx = 1

        for j in 1:ncols
            for i in 1:nrows
                if isnothing(start_date)
                    tick_arg = (xticks = LinearTicks(nxticks),)
                else
                    tick_arg = NamedTuple()
                end
                make_ax(pos) = Axis(pos;
                    xminorticksvisible = true,
                    xminorgridvisible = true,
                    xminorticks = IntervalsBetween(5),
                    yminorticksvisible = true,
                    yminorgridvisible = true,
                    yminorticks = IntervalsBetween(5),
                    tick_arg...
                )
                if selectors
                    plot_box = GridLayout(plot_layout[i, j], 2, 2)
                    well_or_fld, name = split_name(plots[plot_idx])
                    submenu1, l1 = label_menu(plot_box[1, 1], source_keys, "Source", default = well_or_fld)
                    submenu2, l2 = label_menu(plot_box[1, 2], get_quantity_options(well_or_fld), "Type", default = name)
                    # Capture variable in local scope for the functions
                    local_plot_idx = plot_idx
                    on(submenu1.selection) do s
                        update_field_or_well_selection(submenu1, submenu2, s, local_plot_idx)
                    end
                    on(submenu2.selection) do s
                        update_response_selection(submenu1, submenu2, s, local_plot_idx)
                    end
                    subax = make_ax(plot_box[2, 1:2])
                else
                    plot_box = GridLayout(plot_layout[i, j], 1, 1)
                    subax = make_ax(plot_box[1, 1])
                    submenu1 = submenu2 = l1 = l2 = missing
                end
                if hideaxes
                    hidexdecorations!(subax)
                    hideydecorations!(subax)
                elseif i == nrows
                    if isnothing(start_date)
                        subax.xlabel[] = "days"
                    end
                else
                    hidexdecorations!(subax)
                end

                plot_boxes[i, j] = (ax = subax, menu1 = submenu1, menu2 = submenu2, box = plot_box, label1 = l1, label2 = l2)
                plot_idx += 1
            end
        end
        linkxaxes!(map(x -> x.ax, plot_boxes)...)
        update_plots()
    end

    on(row_menu.selection) do _
        update_menu_layout()
    end
    on(col_menu.selection) do _
        update_menu_layout()
    end
    on(unit_menu.selection) do _
        for i in eachindex(plots)
            update_plots(i)
        end
    end
    update_menu_layout()
    if !ismissing(events)
        DataInspector(fig)
    end

    function cycle_well_or_field(inc = 1)
        for pb in plot_boxes
            m1 = pb.menu1
            m2 = pb.menu2
            curr = m1.selection[]
            if curr != "FIELD"
                opts = m1.options[]
                pos = findfirst(isequal(curr), opts)
                pos = mod1(pos + inc, length(opts))
                if pos == 1
                    pos = 2
                end
                m1.selection[] = opts[pos]
                m1.i_selected[] = pos
            end
        end
    end

    function cycle_response(inc = 1)
        for pb in plot_boxes
            m1 = pb.menu1
            m2 = pb.menu2
            curr = m2.selection[]
            opts = m2.options[]
            pos = findfirst(isequal(curr), opts)
            if isnothing(pos)
                pos = 1
            end
            pos = mod1(pos + inc, length(opts))
            m2.selection[] = opts[pos]
            m2.i_selected[] = pos
        end
    end
    on(Makie.events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat
            if event.key == Keyboard.left
                # Previous response
                cycle_response(-1)
            elseif event.key == Keyboard.right
                # Next response
                cycle_response(1)
            elseif event.key == Keyboard.up
                # Previous well/field
                cycle_well_or_field(-1)
            elseif event.key == Keyboard.down
                # Next well/field
                cycle_well_or_field(1)
            end
        end
    end


    return fig
end

function label_menu(dest, options, mlabel::String; kwarg...)
    g = GridLayout(dest)
    l = Label(
        g[1, 1], mlabel,
        justification = :center,
        lineheight = 0.9
    )
    m = Menu(g[1, 2], options = options; kwarg...)
    return (m, l)
end

function plot_events!(ax, events, well_or_fld, start_date, t, v)
    events_for_k = get(events, well_or_fld, missing)
    if ismissing(events_for_k)
        return ax
    end
    ed = collect(keys(events_for_k))
    vals = collect(values(events_for_k))
    if length(ed) == 0
        return ax
    end
    if first(ed) isa DateTime

        t_s = map(x -> Dates.value(x - start_date) / 1000, t)
        event_s = map(x -> Dates.value(x - start_date) / 1000, ed)
    else
        t_s = t
        event_s = ed
    end
    vals_h = get_1d_interpolator(t_s, v).(event_s)
    function event_inspector(plot, index, position)
        return "$(ed[index])\n$(vals[index])"
    end

    scatter!(ax, ed, vals_h; color = :darkred, markersize = 3.0, inspector_label = event_inspector)
    return ax
end
