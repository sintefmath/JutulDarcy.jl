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
        geometry = tpfv_geometry(g)
        centers = geometry.cell_centroids
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
        txt = text!(ax, well_name_for_plot(w, name),
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
            @assert nt == ns "Series $i: Recieved $nt steps, but wells had $ns results."
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
    @assert unit_sys in ("Metric", "SI", "Field")
    t = info.unit_type
    u = 1.0
    if unit_sys == "Metric"
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
    elseif unit_sys == "Field"
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
    elseif unit_sys == "SI"
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

