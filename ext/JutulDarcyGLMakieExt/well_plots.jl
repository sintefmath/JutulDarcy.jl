function JutulDarcy.plot_well!(ax, g, w; color = :darkred, textcolor = nothing, name = nothing, linewidth = 5, top_factor = 0.2, fontsize = 18, geometry = tpfv_geometry(g), kwarg...)
    if isnothing(textcolor)
        textcolor = color
    end
    centers = geometry.cell_centroids
    coord_range(i) = maximum(centers[i, :]) - minimum(centers[i, :])

    if size(centers, 1) == 3
        z = centers[3, :]
    else
        z = [0.0, 1.0]
    end
    bottom = maximum(z)
    top = minimum(z)

    # xrng = coord_range(1)
    # yrng = coord_range(2)

    rng = top - bottom
    s = top + top_factor*rng

    c = well_cells_for_plot(w)
    pts = centers[:, [c[1], c...]]
    if size(pts, 1) == 2
        # 2D grid, add some zeros to make logic work
        pts = vcat(pts, zeros(1, size(pts, 2)))
    end
    pts[3, 1] = s

    l = pts[:, 1]
    text!(well_name_for_plot(w, name),
            position = Tuple([l[1], l[2], -l[3]]),
            space = :data,
            color = textcolor,
            align = (:center, :baseline),
            fontsize = fontsize)
    lines!(ax, vec(pts[1, :]), vec(pts[2, :]), -vec(pts[3, :]), linewidth = linewidth, color = color, kwarg...)
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

# export plot_well_results
function JutulDarcy.plot_well_results(well_data::Dict, arg...; name = "Data", kwarg...)
    JutulDarcy.plot_well_results([well_data], arg...; names = [name], kwarg...)
end

function JutulDarcy.plot_well_results(well_data::Vector, time = nothing; start_date = nothing,
                                                              names =["Dataset $i" for i in 1:length(well_data)], 
                                                              linewidth = 3,
                                                              cmap = nothing, 
                                                              dashwidth = 1,
                                                              new_window = false,
                                                              styles = [:solid, :scatter, :dash, :dashdot, :dot, :dashdotdot],
                                                              resolution = (1600, 900),
                                                              kwarg...)
    # Figure part
    names = Vector{String}(names)
    ndata = length(well_data)
    wd = first(well_data)
    # Selected well
    wells = sort!(collect(keys(wd)))
    nw = length(wells)
    if nw == 0
        return nothing
    end

    is_inj = is_injectors(wd)
    @assert ndata <= length(styles) "Can't plot more datasets than styles provided"
    fig = Figure(resolution = resolution)
    no_time = isnothing(time)
    if no_time
        t_l = "Time-step"
        @assert isnothing(start_date) "start_date does not make sense in the absence of time-steps"
    elseif isnothing(start_date)
        t_l = "Time [days]"
    else
        t_l = "Date"
    end
    ax = Axis(fig[1, 1], xlabel = t_l)

    if isnothing(cmap)
        if nw > 20
            c_key = :turbo
        elseif nw > 10
            c_key = :tab20
        else
            c_key = :Paired_10
        end
        cmap = cgrad(c_key, nw, categorical=true)
    end
    wellstr = [String(x) for x in wells]

    # Type of plot (bhp, rate...)
    responses = collect(keys(wd[first(wells)]))
    respstr = [String(x) for x in responses]
    response_ix = Observable(1)
    is_accum = Observable(false)
    is_abs = Observable(false)
    type_menu = Menu(fig, options = respstr, prompt = respstr[1])

    on(type_menu.selection) do s
        val = findfirst(isequal(s), respstr)
        response_ix[] = val
        autolimits!(ax)
    end
    use_two_cols = nw > 5
    if use_two_cols
        right_block = 2:3
    else
        right_block = 2:2
    end
    fig[2, right_block] = hgrid!(
        type_menu)

    b_xlim = Button(fig, label = "Reset x")
    on(b_xlim.clicks) do n
        reset_limits!(ax; xauto = true, yauto = false)
    end
    b_ylim = Button(fig, label = "Reset y")
    on(b_ylim.clicks) do n
        reset_limits!(ax; xauto = false, yauto = true)
    end
    toggle_abs = Toggle(fig, active = false)
    connect!(is_abs, toggle_abs.active)
    toggle_accum = Toggle(fig, active = false)
    connect!(is_accum, toggle_accum.active)

    buttongrid = GridLayout(tellwidth = false)
    buttongrid[1, 1] = toggle_abs
    buttongrid[1, 2] = Label(fig, "Absolute")
    buttongrid[1, 3] = toggle_accum
    buttongrid[1, 4] = Label(fig, "Cumulative")
    buttongrid[1, 5] = b_xlim
    buttongrid[1, 6] = b_ylim

    # Lay out and do plotting
    fig[2, 1] = buttongrid
    function get_data(time, wix, rix, dataix, use_accum, use_abs)
        tmp = well_data[dataix][wells[wix]][responses[rix]]
        if use_accum && respstr[dataix] != "Bottom hole pressure"
            T = [0.0, time[dataix]...]
            tmp = cumsum(tmp.*diff(T))
        end
        if use_abs
            tmp = abs.(tmp)
        end
        return tmp
    end

    sample = map(x -> get_data([], 1, 1, x, false, false), 1:ndata)
    nsample = map(length, sample)
    if isnothing(time)
        time = map(s -> collect(1:length(s)), sample)
    else
        if eltype(time)<:AbstractFloat# || eltype(time)<:Date
            time = repeat([time], ndata)
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
                    @. T /= (3600*24)
                else
                    T = @. Microsecond(ceil(T*1e6)) + start_date
                end
                push!(newtime, T)
            end
        end
        time = newtime
    end
    lighten(x) = GLMakie.ARGB(x.r, x.g, x.b, 0.2)
    toggles = Vector{Any}([Toggle(fig, active = true, buttoncolor = cmap[i], framecolor_active = lighten(cmap[i])) for i in eachindex(wells)])

    b_inj_on = Button(fig, label = "✔️ I")
    b_inj_off = Button(fig, label = "❌ I")

    labels = Vector{Any}([Label(fig, w) for w in wellstr])

    b_prod_on = Button(fig, label = "✔️ P")
    b_prod_off = Button(fig, label = "❌ P")

    buttongrid[1, 7] = b_inj_on
    buttongrid[1, 8] = b_inj_off
    buttongrid[1, 9] = b_prod_on
    buttongrid[1, 10] = b_prod_off


    function toggle_wells(do_injectors, status)
        for (i, w) in enumerate(wellstr)
            if is_inj[Symbol(w)] == do_injectors
                toggles[i].active[] = status
            end
        end
    end
    on(b_inj_on.clicks) do n
        toggle_wells(true, true)
    end
    on(b_inj_off.clicks) do n
        toggle_wells(true, false)
    end
    on(b_prod_on.clicks) do n
        toggle_wells(false, true)
    end
    on(b_prod_off.clicks) do n
        toggle_wells(false, false)
    end
    tmp = hcat(toggles, labels)
    bgrid = tmp
    N = size(bgrid, 1)

    if use_two_cols
        M = div(N, 2, RoundUp)
        fig[1, 2] = grid!(bgrid[1:M, :], tellheight = false)
        fig[1, 3] = grid!(bgrid[(M+1):N, :], tellheight = false)
    else
        fig[1, 2] = grid!(bgrid, tellheight = false)
    end

    lineh = []
    for dix = 1:ndata
        T = time[dix]
        for i in 1:nw
            d = @lift(get_data(time, i, $response_ix, dix, $is_accum, $is_abs))
            style = styles[dix]
            if style == :scatter
                h = scatter!(ax, T, d, color = cmap[i], linewidth = linewidth, marker = :circle)
            else
                if style == :dash || style == :dashdot
                    lw = dashwidth
                else
                    lw = linewidth
                end
                h = lines!(ax, T, d, linewidth = linewidth, linestyle = style, color = cmap[i])
            end
            t = toggles[i]
            connect!(h.visible, t.active)
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
        fig[1, 1] = Legend(fig, elems, names,
                        tellheight = false,
                        tellwidth = false,
                        margin = (10, 10, 10, 10),
                        halign = :left, valign = :top, orientation = :horizontal
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

