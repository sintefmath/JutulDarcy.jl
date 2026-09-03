"""
    plot_well!(ax, mesh, w; color = :darkred)

Plot a given well that exists in mesh in Axis.
"""
function plot_well!

end

"""
    plot_well_results(wr::WellResults)
    plot_well_results(v::Vector{WellResults})

Launch interactive viewer for well results. Needs GLMakie to be loaded.
"""
function plot_well_results

end


export plot_reservoir_measurables
"""
    plot_reservoir_measurables(case::JutulCase, result::ReservoirSimResult)

Launch interactive viewer for reservoir measurables. Needs GLMakie to be loaded.
"""
function plot_reservoir_measurables

end

export plot_summary
"""
    plot_summary(summary::Dict)
    plot_summary(res::ReservoirSimResult)
    plot_summary([res1, res2, res3]; names = ["Res1", "Res2", "Res3"])

Plot summary results interactively. If multiple results are given, they will be
compared in the same figure.

# Keyword arguments:
- `names`: Names for the different results when multiple results are given.
- `unit_system`: Unit system to use (can be changed in dropdown). Can be
  `:metric`, `:si` or `:field`.
- `linewidth`: Linewidth for the plots.
- `plots` = Vector{Symbol}: Which plots to show by default. For example,
  `["FOPR", "FWPR"]` will field show oil and water production rates as two
  plots. Alternatively, combined plots can be made: `["FOPR,FWPR"]` will show
  both oil and water rates in the same plot. For wells, the name must be
  specified: `["W1:WBHP", "W2:WBHP"]` will show bottom hole pressures for wells
  named W1 and W2.
- `cols`::Int: Number of columns in the layout.
- `selectors`::Bool: Whether to show dropdown selectors for choosing which plots
  to show.
- `extra_field/extra_well`: Additional reservoir measurables or well measurables
  to include in the selection lists, for example to add custom composite plots.
  For example, adding `WBHP,WWIR` to `extra_well` will allow plotting bottom
  hole pressure together with water injection rate for wells.

# Keyboard shortcuts
- `Left Arrow`: Previous response
- `Right Arrow`: Next response
- `Up Arrow`: Previous well (not active for field views)
- `Down Arrow`: Next well (not active for field views)
"""
function plot_summary(arg...; kwarg...)
    Jutul.check_plotting_availability()
    return plot_summary_impl(arg...; kwarg...)
end

function plot_summary_impl

end

function plot_mismatch

end

"""
    plot_reservoir_simulation_result(model::MultiModel, res::ReservoirSimResult; wells = true, reservoir = true)

Plot a reservoir simulation result. If `wells=true` well curves will be shown
interactively. If `reservoir=true` the reservoir quantities will be visualized
in 3D. These options can be combined.
"""
function plot_reservoir_simulation_result(model::MultiModel, res::ReservoirSimResult;
        wells = true,
        reservoir = true,
        fancy = true
    )
    Jutul.check_plotting_availability()
    if reservoir
        fig = plot_reservoir(model, res, fancy = fancy)
    else
        fig = nothing
    end
    if wells
        plot_well_results(res.wells, res.time, new_window = true)
    end
    if reservoir
        display(fig)
    end
    return fig
end

"""
    plot_reservoir(case)
    plot_reservoir(case, states)
    plot_reservoir(case, result)
    plot_reservoir(model, states)

Launch interactive plotter of reservoir + well trajectories in reservoir.
Requires GLMakie to be loaded (using GLMakie). If the keyword `fancy=true`, a
more advanced GUI with more options will be launched that allows for panning and
zooming. If `fancy=false`, a fixed-axis plot will be launched instead. The
keyword `gui=false` can be used to just get a static plot without interactivity
(will use same fixed-axis plot as `fancy=false`).

Displayed units are determined by the `unit_system` keyword, which defaults to
`"metric"`. The source units are assumed to be `"si"` by default. Unit
conversion can be disabled altogether by setting `convert_units=false`.

# Positional arguments
The first entry can be a `model::MultiModel` (from `setup_reservoir_model`), a
`JutulCase` or `reservoir::DataDomain` (from `reservoir_domain`)
"""
function plot_reservoir

end

function plot_reservoir(model, states = missing;
        gui = true,
        fancy = true,
        sens = missing,
        add_secondary = false,
        faults = fancy,
        fault_alpha = 0.5,
        well_fontsize = 18,
        well_linewidth = 3,
        well_color = :darkred,
        zaspect = 1/3,
        aspect = missing,
        well_top_factor_scale = 1.0,
        well_arg = NamedTuple(),
        force_glmakie = true,
        convert_units = true,
        unit_system = "metric",
        source_unit_system = "si",
        unit_lookup = Dict(),
        wells = missing,
        kwarg...
    )
    function maybe_convert_units(x)
        return convert_for_plotting(x,
            units = unit_system,
            source_units = source_unit_system,
            unit_lookup = unit_lookup,
            convert_units = convert_units
        )
    end
    if states isa AbstractDict
        states = [states]
    end
    if model isa DataDomain
        data_domain = model
        model = missing
    else
        data_domain = reservoir_domain(model)
    end

    if !ismissing(model)
        if !ismissing(states) && add_secondary
            rmodel = reservoir_model(model)
            states = [Jutul.evaluate_all_secondary_variables(rmodel, s) for s in states]
        end
    end
    states = maybe_convert_units(states)
    sens = maybe_convert_units(sens)
    Jutul.check_plotting_availability()
    if force_glmakie
        @assert Jutul.plotting_check_interactive(warn = true) "Function requires interactive plotting. Set force_glmakie = false to override."
    end
    cell_centroids = data_domain[:cell_centroids]
    if ismissing(aspect)
        x = cell_centroids[1, :]
        y = cell_centroids[2, :]
        xrng = maximum(x) - minimum(x)
        yrng = maximum(y) - minimum(y)
        aspect = (1.0, max(yrng/xrng, 0.001), zaspect)
    end
    if haskey(data_domain, :boundary_centroids)
        bc = data_domain[:boundary_centroids]
        if size(bc, 1) == 3
            zb = data_domain[:boundary_centroids][3, :]
            filter!(isfinite, zb)
            if length(zb) > 1
                bounds_z = (minimum(zb), maximum(zb))
            else
                bounds_z = missing
            end
        else
            bounds_z = missing
        end
    else
        bounds_z = missing
    end
    g = physical_representation(data_domain)
    static_data = maybe_convert_units(data_domain)

    wtoggle = ftoggle = missing
    if gui
        if fancy
            if !ismissing(aspect)
                aspect = 1.0 ./ aspect
            end
            # In case it is not yet supported in Jutul...
            if ismissing(sens)
                extra_arg = tuple()
            else
                extra_arg = (sens = sens, )
            end
            s = plot_explorer(g;
                dynamic = states,
                static = static_data,
                zreversed = true,
                aspect = aspect,
                extra_arg...,
                kwarg...
            )
            ax = s.lscene
            fig = s.fig
            wtoggle = s.add_toggle("Wells", true)
            ftoggle = s.add_toggle("Faults", faults)
            retval = s
        else
            if ismissing(states)
                arg = tuple()
            else
                merged_states = map(states) do s
                    o = OrderedDict()
                    for (k, v) in pairs(s)
                        o[k] = v
                    end
                    for (k, v) in pairs(static_data)
                        o[k] = v
                    end
                    return o
                end
                arg = (merged_states, )
            end
            fig = plot_interactive(data_domain, arg...; z_is_depth = true, aspect = aspect, kwarg...)
            ax = fig.current_axis[]
            retval = fig
        end
    else
        fig, ax, _ = plot_cell_data(g, arg...; z_is_depth = true, kwarg...)
        retval = fig
    end
    if ismissing(wells)
        wells = Dict{Symbol, Any}()
        if model isa MultiModel
            for (k, m) in pairs(model.models)
                w = physical_representation(m.data_domain)
                if w isa WellDomain
                    wells[k] = w
                end
            end
        end
    elseif wells isa AbstractVector
        ws = Dict{Symbol, Any}()
        for w in wells
            w = physical_representation(w)
            ws[w.name] = w
        end
        wells = ws
    end

    i = 1
    n = length(wells)
    well_plts = []
    for (k, w) in pairs(wells)
        tf = 0.2 + 0.1*(i/n)
        if well_color isa AbstractDict
            well_color_k = get(well_color, k, :darkred)
        else
            well_color_k = well_color
        end
        wp = plot_well!(ax.scene, g, w;
            fontsize = well_fontsize,
            top_factor = well_top_factor_scale*tf,
            bounds_z = bounds_z,
            color = well_color_k,
            linewidth = well_linewidth,
            cell_centroids = cell_centroids,
            extra_out = true,
            toggle = wtoggle,
            well_arg...
        )
        push!(well_plts, wp)
        i += 1
    end
    if faults
        plot_faults!(ax, g; domain = data_domain, toggle = ftoggle, alpha = fault_alpha)
    end
    return retval
end

function plot_reservoir(model::Union{MultiModel, SimulationModel,}, result::ReservoirSimResult; kwarg...)
    return plot_reservoir(model, result.states; kwarg...)
end

function plot_reservoir(case::JutulCase, states = missing; state0 = true, kwarg...)
    if states isa ReservoirSimResult
        states = states.states
    end
    if states isa AbstractDict
        states = [states]
    end
    if state0
        s0 = case.state0
        m = case.model
        if haskey(s0, :Reservoir)
            s0 = s0[:Reservoir]
            m = m[:Reservoir]
        end
        s0 = Jutul.evaluate_all_secondary_variables(m, s0)
        if ismissing(states) || length(states) == 0
            states = [s0]
        else
            s = states[1]
            s0_new = typeof(s)()
            for (k, v) in pairs(s)
                s0_new[k] = v
            end
            states = [s0_new; states]
        end
    end
    return plot_reservoir(case.model, states; kwarg...)
end

export plot_faults!
function plot_faults!(ax, domain::DataDomain; kwarg...)
    return plot_faults!(ax, physical_representation(domain); domain = domain, kwarg...)
end

"""
    fig = JutulDarcy.plot_co2_inventory(t, inventory, plot_type = :stack)

Plots the CO2 inventory over time or steps, with options for stacked or line
plots. `inventory` is the output from `co2_inventory` while `t` can either be
omitted, be a list of reporting time in seconds or a index list of steps where
the solution is given.

# Arguments
- `t`: A vector representing time or steps. If `t` is of type `Float64`, it is
  assumed to represent time in seconds and will be converted to years.
- `inventory`: A vector of dictionaries, where each dictionary contains CO2 mass
  data for different categories (e.g., `:dissolved`, `:mobile`, `:residual`,
  etc.).
- `plot_type`: (Optional) A symbol specifying the type of plot. Can be `:stack`
  for stacked plots or `:lines` for line plots. Default is `:stack`.

# Notes

This function is only available if Makie is loaded (through for example GLMakie
or CairoMakie)
"""
function plot_co2_inventory

end

"""
    plot_well_vs_meters_drilled!(ax, well_model, values; label=nothing, kwargs...)

Plot well data vs meters drilled on the given axis. If the well has sections,
each section is plotted with a distinct color. Needs GLMakie to be loaded.
"""
function plot_well_vs_meters_drilled!

end

"""
    plot_well_vs_meters_drilled(well_model, values; figure_kwargs=NamedTuple(), axis_kwargs=NamedTuple(), kwargs...)

Create a new figure and plot well data vs meters drilled. Needs GLMakie to be loaded.
"""
function plot_well_vs_meters_drilled

end

"""
    plot_well_states_interactive(well_model, states; time=missing, names=missing, resolution=(1200, 800), kwargs...)
    plot_well_states_interactive(well_name::Symbol, model, states; kwargs...)

Interactive GUI for plotting multisegment well data vs meters drilled. Needs GLMakie to be loaded.
"""
function plot_well_states_interactive

end
