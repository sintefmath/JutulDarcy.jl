using Dates

function JutulDarcy.plot_well!(ax, g, w;
        color = :darkred,
        textcolor = nothing,
        name = nothing,
        linewidth = 3,
        top_factor = 0.2,
        fontsize = 18,
        glowwidth = 5.0,
        markersize = 10.0,
        transparency = true,
        bounds_z = missing,
        cell_centroids = missing,
        kwarg...
    )
    if isnothing(textcolor)
        textcolor = color
    end
    c = well_cells_for_plot(w)
    if ismissing(cell_centroids)
        if g isa DataDomain
            centers = g[:cell_centroids]
            if ismissing(bounds_z) && haskey(g, :boundary_centroids)
                fc = g[:boundary_centroids]
                if size(fc, 1) == 3
                    bounds_z = extrema(view(fc, 3, :))
                end
            end
        else
            geometry = tpfv_geometry(g)
            centers = geometry.cell_centroids
        end
    else
        centers = cell_centroids
    end
    if size(centers, 1) == 3
        z = view(centers, 3, :)
    else
        z = [0.0, 1.0]
    end
    if ismissing(bounds_z)
        bottom = maximum(z)
        top = minimum(z)
    else
        top, bottom = bounds_z
    end

    rng = top - bottom
    s = top + top_factor*rng

    pts = centers[:, [c[1], c...]]
    if size(pts, 1) == 2
        # 2D grid, add some zeros to make logic work
        pts = vcat(pts, zeros(1, size(pts, 2)))
    end
    pts[3, 1] = s

    l = pts[:, 1]
    if fontsize > 0
        wd = physical_representation(w)
        txt = text!(ax, well_name_for_plot(wd, name),
            position = Tuple([l[1], l[2], l[3]]),
            space = :data,
            transparency = transparency,
            color = textcolor,
            align = (:center, :baseline),
            glowcolor = (:white, 1.0),
            glowwidth = glowwidth,
            fontsize = fontsize
        )
    end
    x = pts[1, :]
    y = pts[2, :]
    z = pts[3, :]

    npts = length(x)
    spts = fill(markersize, npts)
    spts[1] = 0
    scatterlines!(ax, x, y, z,
        transparency = transparency,
        linewidth = linewidth,
        color = color,
        markersize = spts,
        alpha = 0.9,
        kwarg...
    )
end

well_name_for_plot(w::Dict, ::Nothing) = w["name"]
well_name_for_plot(w, s::String) = s
well_name_for_plot(w, s::Nothing) = String(w.name)

function well_cells_for_plot(w::Dict)
    wc = w["cells"]
    if !isa(wc, AbstractArray)
        wc = [wc]
    end
    return vec(Int64.(wc))
end

function well_cells_for_plot(w::DataDomain)
    return w |> physical_representation |> well_cells_for_plot
end

function well_cells_for_plot(w)
    return w.perforations.reservoir
end

function JutulDarcy.plot_well_results(r::JutulDarcy.ReservoirSimResult; kwarg...)
    return JutulDarcy.plot_well_results(r.wells; kwarg...)
end

function JutulDarcy.plot_well_results(well_data, arg...; name = "Data", kwarg...)
    JutulDarcy.plot_well_results([well_data], arg...; names = [name], kwarg...)
end

function JutulDarcy.plot_well_results(well_data::Vector, time = missing;
        start_date = first(well_data).start_date,
        names = ["Dataset $i" for i in 1:length(well_data)],
        linewidth = 1.5,
        cmap = nothing,
        dashwidth = 1,
        field = missing,
        accumulated = false,
        unit_sys = "Metric",
        new_window = false,
        styles = [:solid, :dash, :scatter, :dashdot, :dot, :dashdotdot],
        resolution = (1600, 900),
        kwarg...
    )
    LEFT_COLUMN_WIDTH = 180
    LEFT_COLUMN_ITEM_WIDTH = LEFT_COLUMN_WIDTH - 10
    response_ix = Observable(1)
    is_accum = Observable(accumulated)
    is_abs = Observable(false)
    is_line = Observable(true)
    unit_sys = Observable(unit_sys)

    # Figure part
    names = Vector{String}(names)
    ndata = length(well_data)
    wd = first(well_data).wells
    # Selected well
    wells = sort!(collect(keys(wd)))
    nw = length(wells)
    if nw == 0
        return nothing
    end
    # Type of plot (bhp, rate...)
    responses = collect(keys(wd[first(wells)]))
    setdiff!(responses, [:control])
    respstr = [String(x) for x in responses]

    is_inj = is_injectors(wd)
    @assert ndata <= length(styles) "Can't plot more datasets than styles provided"
    fig = Figure(size = resolution)
    no_time = isnothing(time)
    if no_time
        t_l = "Time-step"
        @assert isnothing(start_date) "start_date does not make sense in the absence of time-steps"
    elseif isnothing(start_date)
        t_l = "days"
    else
        t_l = "Date"
    end
    function response_label_to_unit(s, is_accumulated)
        info = JutulDarcy.well_target_information(Symbol(s))
        if ismissing(info)
            return ""
        else
            lbl = info.unit_label
            if is_accumulated[]
                lbl = replace(lbl, "/s" => "")
            end
            return lbl
        end
    end
    function response_label_to_descr(s)
        info = JutulDarcy.well_target_information(Symbol(s))
        if ismissing(info)
            return "$s"
        else
            return "$(info.description)"
        end
    end
    if ismissing(field) || !(field in responses)
        field = first(responses)
    end
    y_l = Observable(response_label_to_unit(field, is_accum))
    title_l = Observable(response_label_to_descr(field))

    ax = Axis(fig[1, 2], xlabel = t_l, ylabel = y_l, title = title_l)

    if isnothing(cmap)
        if nw > 20
            c_key = :turbo
        elseif nw > 10
            c_key = :tab20
        else
            c_key = :Paired_10
        end
        cmap = cgrad(c_key, max(nw, 2), categorical=true)
    end
    wellstr = [String(x) for x in wells]

    type_menu = Menu(fig, options = respstr, prompt = respstr[1], width = LEFT_COLUMN_ITEM_WIDTH)

    on(type_menu.selection) do s
        val = findfirst(isequal(s), respstr)
        y_l[] = response_label_to_unit(s, is_accum)
        title_l[] = response_label_to_descr(s)
        response_ix[] = val
        autolimits!(ax)
    end

    unit_menu = Menu(fig, options = ["Metric", "SI", "Field"], prompt = unit_sys[], width = LEFT_COLUMN_ITEM_WIDTH)

    on(unit_menu.selection) do s
        unit_sys[] = s
        notify(response_ix)
        autolimits!(ax)
    end

    b_xlim = Button(fig, label = "Reset x axis", width = LEFT_COLUMN_ITEM_WIDTH)
    on(b_xlim.clicks) do n
        reset_limits!(ax; xauto = true, yauto = false)
    end
    b_ylim = Button(fig, label = "Reset y axis", width = LEFT_COLUMN_ITEM_WIDTH)
    on(b_ylim.clicks) do n
        reset_limits!(ax; xauto = false, yauto = true)
    end
    buttongrid = GridLayout(width = LEFT_COLUMN_WIDTH, tellheight = false, valign = :top)
    button_ix = 1

    buttongrid[button_ix, 1:2] = Label(
        fig, "Plot selection",
        font = :bold
    )
    # Select well response
    button_ix += 1
    buttongrid[button_ix, 1:2] = Label(
        fig, "Well response"
    )
    button_ix += 1
    buttongrid[button_ix, 1:2] = type_menu
    button_ix += 1
    # Select unit system
    buttongrid[button_ix, 1:2] = Label(
        fig, "Unit system"
    )
    button_ix += 1
    buttongrid[button_ix, 1:2] = unit_menu
    button_ix += 1
    # Reset limits
    buttongrid[button_ix, 1:2] = Label(
        fig, "Reset limits"
    )
    button_ix += 1
    buttongrid[button_ix, 1:2] = b_xlim
    button_ix += 1
    buttongrid[button_ix, 1:2] = b_ylim
    button_ix += 1

    toggle_abs = Checkbox(fig, checked = true)
    connect!(is_abs, toggle_abs.checked)
    buttongrid[button_ix, 1] = toggle_abs
    buttongrid[button_ix, 2] = Label(fig, "Absolute value", halign = :left)
    button_ix += 1

    toggle_accum = Checkbox(fig, checked = is_accum[])
    connect!(is_accum, toggle_accum.checked)
    buttongrid[button_ix, 1] = toggle_accum
    buttongrid[button_ix, 2] = Label(fig, "Cumulative sum", halign = :left)
    button_ix += 1

    # toggle_line = Checkbox(fig, checked = true)
    # connect!(is_line, toggle_line.checked)
    # buttongrid[button_ix, 1] = toggle_line
    # buttongrid[button_ix, 2] = Label(fig, "Line", halign = :left)
    # button_ix += 1

    # Lay out and do plotting
    fig[1, 3] = buttongrid
    function get_data(time, wix, rix, dataix, use_accum, use_abs)
        qoi_val = well_data[dataix][wells[wix]][responses[rix]]
        factor = 1.0
        s = responses[response_ix[]]
        lbl = response_label_to_unit(s, use_accum)
        info = JutulDarcy.well_target_information(Symbol(s))

        if ismissing(info)
            is_rate = false
        else
            is_rate = endswith(info.unit_label, "/s")
            lbl = replace(info.unit_label, "/s" => "")
            @assert info.is_rate == is_rate
            f_u, lbl = well_unit_conversion(unit_sys[], lbl, info)
            if f_u isa Symbol
                qoi_val = convert_from_si.(qoi_val, f_u)
            else
                factor *= f_u
            end
        end
        if is_rate
            if !use_accum
                lbl *= "/day"
            end
            factor = factor.*si_unit(:day)
        end
        y_l[] = lbl
        qoi_val = factor.*qoi_val
        if use_accum && is_rate
            T = [0.0, time[dataix]...]
            qoi_val = cumsum(qoi_val.*diff(T))
        end
        if use_abs
            qoi_val = abs.(qoi_val)
        end
        return qoi_val
    end

    sample = map(x -> get_data([], 1, 1, x, false, false), 1:ndata)
    nsample = map(length, sample)
    if isnothing(time)
        time = map(s -> collect(1:length(s)), sample)
    else
        if eltype(time)<:AbstractFloat# || eltype(time)<:Date
            time = repeat([time], ndata)
        else
            time = map(x -> x.time, well_data)
        end
    end
    if !no_time
        newtime = []
        for i = 1:ndata
            T = copy(time[i])
            nt = length(T)
            ns = nsample[i]
            @assert nt == ns "Series $i: Received $nt steps, but wells had $ns results."
            if eltype(T)<:AbstractFloat
                # Scale to days
                if isnothing(start_date)
                    @. T /= si_unit(:day)
                else
                    T = @. Microsecond(ceil(T*1e6)) + start_date
                end
                push!(newtime, T)
            end
        end
        time = newtime
    end
    lighten(x) = GLMakie.ARGB(x.r, x.g, x.b, 0.2)
    toggles = [
            Checkbox(fig,
            checked = true,
            checkboxstrokecolor_checked = lighten(cmap[i]),
            checkboxstrokecolor_unchecked = cmap[i],
            checkboxcolor_checked = cmap[i]
            ) for i in eachindex(wells)
        ]
    # toggles = Vector{Any}([Toggle(fig, active = true, buttoncolor = cmap[i], framecolor_active = lighten(cmap[i])) for i in eachindex(wells)])


    labels = Vector{Any}([Label(fig, w) for w in wellstr])

    tmp = hcat(toggles, labels)
    bgrid = tmp
    N = size(bgrid, 1)

    use_two_cols = nw > 20
    leg_layout = fig[1, 1] = GridLayout()
    wsel_label = Label(
        fig, "Well selection",
        font = :bold
    )
    if use_two_cols
        M = div(N, 2, RoundUp)
        leg_layout[2, 1] = grid!(bgrid[1:M, :], tellheight = false)
        leg_layout[2, 2] = grid!(bgrid[(M+1):N, :], tellheight = false)
    else
        leg_layout[2, 1] = grid!(bgrid, tellheight = false)
    end
    leg_layout[1, :] = wsel_label

    b_prod_on = Button(fig, label = "Producers")
    b_inj_on = Button(fig, label = "Injectors")
    b_all_on = Button(fig, label = "All")

    well_toggle_group = GridLayout(tellwidth = true)
    well_toggle_group[2, 1] = b_prod_on
    well_toggle_group[2, 2] = b_inj_on
    well_toggle_group[2, 3] = b_all_on
    well_toggle_group[1, :] = Label(fig, "Toggle wells")

    leg_layout[3, :] = well_toggle_group

    function toggle_wells(wtype)
        do_injectors = wtype == :injectors
        do_all = wtype == :all
        num_active = 0
        crit(w) = is_inj[Symbol(w)] == do_injectors || do_all
        for (i, w) in enumerate(wellstr)
            if crit(w)
                num_active += toggles[i].checked[]
            end
        end
        status = num_active == 0
        for (i, w) in enumerate(wellstr)
            if crit(w)
                toggles[i].checked[] = status
            end
        end
    end
    on(b_inj_on.clicks) do n
        toggle_wells(:injectors)
    end
    on(b_prod_on.clicks) do n
        toggle_wells(:producers)
    end
    on(b_all_on.clicks) do n
        toggle_wells(:all)
    end

    lineh = []
    for dix = 1:ndata
        T = time[dix]
        for i in 1:nw
            d = @lift(get_data(time, i, $response_ix, dix, $is_accum, $is_abs))
            style = styles[dix]
            if style == :scatter
                h = scatter!(ax, T, d, color = cmap[i], marker = :circle)
            else
                if style == :dash || style == :dashdot
                    lw = dashwidth
                else
                    lw = linewidth
                end
                lw = @lift $is_line*linewidth

                h = lines!(ax, T, d,
                    linewidth = lw,
                    linestyle = style,
                    color = cmap[i]
                )
            end
            t = toggles[i]
            connect!(h.visible, t.checked)
            push!(lineh, h)
        end
    end
    x = (tmp, tmp2) -> autolimits!(ax)
    @lift(x($is_abs, $is_accum))
    if ndata > 1
        elems = []
        for i = 1:ndata
            style = styles[i]
            if style == :scatter
                el = MarkerElement(color = :black, linewidth = linewidth, marker = :circle)
            else
                if style == :dash || style == :dashdot
                    lw = dashwidth
                else
                    lw = linewidth
                end
                el = LineElement(color = :black, linestyle = style, linewidth = lw)
            end
            push!(elems, el)
        end
        fig[1, 2] = Legend(fig, elems, names,
            tellheight = false,
            tellwidth = false,
            margin = (10, 10, 10, 10),
            halign = :left,
            valign = :top,
            orientation = :horizontal
        )
    end
    if new_window
        display(GLMakie.Screen(), fig)
    end
    return fig
end

function is_injectors(well_data)
    lkey = Symbol("Surface total rate")
    skey = :rate
    D = Dict{Symbol, Bool}()
    for (k, v) in well_data
        if haskey(v, lkey)
            val = v[lkey]
        else
            val = v[skey]
        end
        D[k] = sum(val) > 0
    end
    return D
end

function well_unit_conversion(unit_sys, lbl, info)
    unit_sys = lowercase(unit_sys)
    @assert unit_sys in ("metric", "si", "field")
    t = info.unit_type
    u = 1.0
    if unit_sys == "metric"
        if t == :gas_volume_surface
            lbl = "m³"
        elseif t == :gas_volume_reservoir
            lbl = "m³"
        elseif t == :liquid_volume_surface
            lbl = "m³"
        elseif t == :liquid_volume_reservoir
            lbl = "m³"
        elseif t == :pressure
            u = :bar
            lbl = "bar"
        elseif t == :absolute_temperature
            u = :Celsius
            lbl = "°C"
        elseif t == :relative_temperature
            u = :Kelvin
            lbl = "°K"
        end
    elseif unit_sys == "field"
        if t == :gas_volume_surface
            u = si_unit(:kilo)*si_unit(:feet)^3
            lbl = "MScf"
        elseif t == :gas_volume_reservoir
            u = :stb
            lbl = "bbl"
        elseif t == :liquid_volume_surface
            u = :stb
            lbl = "bbl"
        elseif t == :liquid_volume_reservoir
            u = :stb
            lbl = "bbl"
        elseif t == :pressure
            u = :psi
            lbl = "psi"
        elseif t == :absolute_temperature
            u = :Rankine
            lbl = "°R"
        elseif t == :relative_temperature
            u = :Fahrenheit
            lbl = "°F"
        elseif t == :mass
            u = :pound
            lbl = "pound"
        end
    elseif unit_sys == "si"
        if t == :gas_volume_surface
            lbl = "m³"
        elseif t == :gas_volume_reservoir
            lbl = "m³"
        elseif t == :liquid_volume_surface
            lbl = "m³"
        elseif t == :liquid_volume_reservoir
            lbl = "m³"
        elseif t == :pressure
            lbl = "Pa"
        elseif t == :absolute_temperature
            lbl = "°K"
        elseif t == :relative_temperature
            lbl = "°C"
        elseif t == :mass
            lbl = "kg"
        end
    end
    return (u, lbl)
end

function JutulDarcy.plot_reservoir_measurables(arg...;
        type = :field,
        left = missing,
        right = "none",
        lcolor = Makie.wong_colors()[1],
        rcolor = Makie.wong_colors()[6],
        accumulated = false,
        left_accumulated = accumulated,
        right_accumulated = accumulated,
        unit_system = "Metric",
        kwarg...
    )
    if length(arg) == 1
        fieldvals = only(arg)
    else
        fieldvals = JutulDarcy.reservoir_measurables(arg..., type = type)
    end
    fig = Figure(size = (1200, 800))

    t = fieldvals[:time]/si_unit(:day)
    dt = diff([0.0, t...])
    n = length(t)
    function get_field_data(ax, k, accumulated, usel)
        if k == "none"
            out = fill(NaN, n)
            u = ""
        else
            data = fieldvals[to_key(k)]
            out = data.values
            factor, u = well_unit_conversion(usel, "", data)
            if factor isa Symbol
                out = convert_from_si.(out, factor)
            else
                factor *= factor
            end
            if data.is_rate
                out = out.*si_unit(:day)
                if accumulated
                    out = cumsum(out .* dt)
                else
                    u *= "/day"
                end
            end
        end
        ax.ylabel[] = u
        return out
    end

    mkeys = String[]
    left_found = right_found = false
    for (k, v) in pairs(fieldvals)
        if k == :time
            continue
        end
        next = "$k - $(fieldvals[k].legend)"
        push!(mkeys, next)
        if !ismissing(left) && "$k" == "$left"
            left = next
            left_found = true
        end
        if !ismissing(right) && "$k" == "$right"
            right = next
            right_found = true
        end
    end
    sort!(mkeys)
    push!(mkeys, "none")
    if ismissing(right) || !right_found
        right = mkeys[1]
    end
    if ismissing(left) || !left_found
        left = mkeys[1]
    end
    bg = GridLayout(tellheight = true)
    fig[1, 1] = bg
    # Plotting axis
    ax1 = Axis(fig[2, 1],
        yticklabelcolor = lcolor,
        xlabel = "days",
        ylabelcolor = lcolor
    )
    ax2 = Axis(fig[2, 1],
        yticklabelcolor = rcolor,
        yaxisposition = :right,
        ylabelcolor = rcolor
    )
    hidespines!(ax2)
    hidexdecorations!(ax2)
    ax2.xgridvisible = false
    ax2.ygridvisible = false
    deactivate_interaction!(ax2, :rectanglezoom)

    function to_key(x::String)
        xs, = split(x, " - ")
        return Symbol(xs)
    end

    function selection_function(sel, s, ax)
        sel[] = s
        if s == "none"
            ylims!(ax, (0.0, 1.0))
        end
        autolimits!(ax)
    end

    function setup_menu(ax, default, pos, color, is_accumulated)
        selection = Observable(default)
        lmenu = Menu(fig, default = selection[], options = mkeys, tellwidth = false)
        if default == "none"
            init = ""
        else
            init = fieldvals[to_key(selection[])].legend
        end
        l_accum = Checkbox(fig, checked = is_accumulated)
        is_accum = Observable(is_accumulated)
        connect!(is_accum, l_accum.checked)

        on(lmenu.selection) do s
            selection_function(selection, s, ax)
        end
        bg[1, pos] = lmenu
        lgroup = GridLayout()
        unit_menu = Menu(fig,
            options = ["Metric", "SI", "Field"],
            prompt = unit_system,
        )
        lgroup[1, 1] = unit_menu
        reset_left = Button(fig, label = "Reset axis")
        on(reset_left.clicks) do _
            autolimits!(ax)
        end

        lgroup[1, 2] = reset_left
        lgroup[1, 3] = l_accum
        lgroup[1, 4] = Label(fig, "Accumulated")
        bg[2, pos] = lgroup

        # Actual plotting
        usel = unit_menu.selection
        d = @lift(get_field_data(ax, $selection, $is_accum, $usel))
        lines!(ax, t, d; color = color, kwarg...)
    end

    # Left side menu
    setup_menu(ax1, left, 1, lcolor, left_accumulated)

    # Right side menu
    setup_menu(ax2, right, 2, rcolor, right_accumulated)

    return fig
end

function JutulDarcy.plot_summary(v::Union{AbstractVector, Tuple}; kwarg...)
    return JutulDarcy.plot_summary(v...; kwarg...)
end

function JutulDarcy.plot_summary(arg...;
        names = ["Dataset $i" for i in 1:length(arg)],
        unit_system = "Metric",
        linewidth = 2.0,
        markersize = 0.0,
        plots = ["FIELD:FPR"],
        extra_field = String[],
        extra_well = String[],
        linecolor = :black,
        selectors = true,
        cols = 1,
        alpha = 1.0,
        rows = ceil(length(plots)/cols) |> Int,
        colormap = missing,
        kwarg...
    )
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
            out = start_date .+ @. Microsecond(ceil(t*1e6))
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
        if !ismissing(info)
            JutulDarcy.GeoEnergyIO.InputParser.swap_unit_system!(v, systems, info.unit_type)
        end
        return v
    end

    # Set up time labels
    start_dates = map(s -> s["TIME"].start_date, summaries)
    unique!(filter!(x -> !isnothing(x) && !ismissing(x), start_dates))
    if length(start_dates) == 0
        start_date = nothing
    else
        if length(start_dates) > 1
            println("Note: Multiple distinct start dates ($start_dates) found in summaries, using first entry.")
        end
        start_date = first(start_dates)
    end
    # Max time - in days
    t_max = maximum(x -> maximum(time_data(x)), eachindex(summaries))
    ticks_days = collect(range(0.0, ceil(t_max), length = 10))
    if !isnothing(start_date)
    end

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
    markersize_menu, = label_menu(top_menu_grid[1, 4], markeroptions, "Marker size", default = string(markersize))
    unit_menu, = label_menu(top_menu_grid[1, 5], ["Metric", "SI", "Field"], "Unit system"; default = unit_system)

    plot_layout = nothing
    plot_boxes = []
    legends = Dict{Int, Legend}()

    function plot_string(source::String, type::String)
        return "$source:$type"
    end

    function plot_string(source::String, ::Nothing)
        return plot_string(source, "NONE")
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
                    if !ismissing(info) && name != "NONE"
                        # ax.title[] = "$(name) $(info.legend)"
                        _, u = well_unit_conversion(units, "", info)
                        if info.is_rate
                            v = v.*si_unit(:day)
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
                    ms_sel = markersize_menu.selection
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
                    scatterlines!(ax, t, v;
                        linewidth = lw_sel,
                        markersize = ms_sel,
                        label = lbl,
                        alpha = alpha,
                        transparency = alpha < 1.0,
                        arg...
                    )
                end
                ax.ylabel[] = ystr
                ax.title[] = tstr
            end
            if nlines > 1
                l = axislegend(ax, position = :lt)
                if haskey(legends, i)
                    delete!(legends[i])
                end
                legends[i] = l
            end
        end
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
                if selectors
                    plot_box = GridLayout(plot_layout[i, j], 2, 2)
                    well_or_fld, name = split_name(plots[plot_idx])
                    submenu1, l1 = label_menu(plot_box[1, 1], source_keys, "Source", default = well_or_fld)
                    submenu2, l2 = label_menu(plot_box[1, 2], get_quantity_options(well_or_fld), "Type", default = name)
                    # Capture variable in local scope for the functions
                    local_plot_idx = plot_idx
                    on(submenu1.selection) do s
                        submenu2.options[] = get_quantity_options(s)
                        plots[local_plot_idx] = plot_string(s, submenu2.selection[])
                        update_plots(local_plot_idx)
                    end
                    on(submenu2.selection) do s
                        plots[local_plot_idx] = plot_string(submenu1.selection[], s)
                        update_plots(local_plot_idx)
                    end
                    subax = Axis(plot_box[2, 1:2])
                else
                    plot_box = GridLayout(plot_layout[i, j], 1, 1)
                    subax = Axis(plot_box[1, 1])
                    submenu1 = submenu2 = l1 = l2 = missing
                end
                if i == nrows
                    if isnothing(start_date) || ncols > 2
                        # subax.xticks[] = ticks_days
                    else
                        dates = map(d -> start_date + Day(round(Int, d)), ticks_days)
                        ticks_dates = map(d -> Dates.format(d, "yyyy u d"), dates)
                        # subax.xticks[] = (ticks_days, ticks_dates)
                    end
                    subax.xlabel[] = "days"
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

function JutulDarcy.plot_mismatch(obj, res::ReservoirSimResult)
    _, w = JutulDarcy.HistoryMatching.evaluate_match(obj, res, log = true)
    return JutulDarcy.plot_mismatch(w)
end

function JutulDarcy.plot_mismatch(w)
    (; impact, targets, wells) = JutulDarcy.HistoryMatching.compute_well_target_contribution_matrix(w)
    nw, ntargets = size(impact)

    fig = Figure(size = (1400, 900))
    # Barplot over well contributions, banked by control type
    ax1 = Axis(fig[2, 1], title = "Contribution by well (total)")
    ax2 = Axis(fig[2, 2], title = "Contribution by target (total)")
    ax3 = Axis(fig[3, 1:2], title = "Contribution by target (per step)")
    ax4 = Axis(fig[4, 1:2], title = "Contribution by well (per step)")

    colors = Makie.wong_colors()

    well_idx = repeat((1:nw), 1, ntargets)
    target_idx = repeat((1:ntargets)', nw, 1)
    tvec = vec(target_idx)
    wvec = vec(well_idx)
    barplot!(ax1, wvec, vec(impact), stack = wvec, color = colors[tvec])

    labels = map(String, targets)
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    Legend(fig[1, 1:2], elements, labels, orientation  = :horizontal, tellheight = true)
    ax1.xticks[] = (1:nw, map(string, wells))
    ax1.xticklabelrotation[] = -pi/4

    # Barplot over control types
    target_weight = sum(impact, dims = 1)
    barplot!(ax2, 1:ntargets, vec(target_weight), color = colors[1:ntargets])
    ax2.xticks[] = (1:ntargets, map(string, targets))

    # Magnitude of control by step
    maxstep = 0
    for well in wells
        for wl in w[well]
            maxstep = max(maxstep, wl.stop)
        end
    end
    target_step_vals = zeros(ntargets, maxstep)
    well_step_vals = zeros(nw, maxstep)
    for step in 1:maxstep
        (; impact, targets, wells) = JutulDarcy.HistoryMatching.compute_well_target_contribution_matrix(w, step = step)
        target_step_vals[:, step] = sum(impact, dims = 1)'
        well_step_vals[:, step] = sum(impact, dims = 2)
    end
    for i in 1:ntargets
        lines!(ax3, 1:maxstep, target_step_vals[i, :], color = colors[i], label = string(targets[i]))
    end
    axislegend(ax3)
    # Magnitude of well by step
    plts = []
    for i in 1:nw
        plt = lines!(ax4, 1:maxstep, well_step_vals[i, :], label = string(wells[i]))
        push!(plts, plt)
    end
    hidexdecorations!(ax3)
    Legend(fig[5, 1:2], plts, map(string, wells), orientation  = :horizontal, tellheight = true, nbanks = 3)

    fig
end
