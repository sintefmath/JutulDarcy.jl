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
"""
function plot_summary

end

function plot_mismatch

end

"""
    plot_reservoir_simulation_result(model::MultiModel, res::ReservoirSimResult; wells = true, reservoir = true)

Plot a reservoir simulation result. If `wells=true` well curves will be shown
interactively. If `reservoir=true` the reservoir quantities will be visualized
in 3D. These options can be combined.
"""
function plot_reservoir_simulation_result(model::MultiModel, res::ReservoirSimResult; wells = true, reservoir = true)
    Jutul.check_plotting_availability()
    if reservoir
        rmodel = reservoir_model(model)
        fig = plot_interactive(rmodel, res.states)
        g = physical_representation(rmodel.data_domain)
        ax = fig.current_axis[]
        for (k, m) in pairs(model.models)
            w = physical_representation(m.data_domain)
            if w isa WellDomain
                plot_well!(ax, g, w)
            end
        end
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
    plot_reservoir(model, states=missing; well_fontsize = 18, well_linewidth = 3, kwarg...)

Launch interactive plotter of reservoir + well trajectories in reservoir. Requires GLMakie.
"""
function plot_reservoir(model, arg...;
        gui = true,
        well_fontsize = 18,
        well_linewidth = 3,
        well_color = :darkred,
        zaspect = 1/3,
        aspect = missing,
        well_top_factor_scale = 1.0,
        well_arg = NamedTuple(),
        force_glmakie = true,
        wells = missing,
        kwarg...
    )
    Jutul.check_plotting_availability()
    if force_glmakie
        @assert Jutul.plotting_check_interactive(warn = true) "Function requires interactive plotting. Set force_glmakie = false to override."
    end
    if model isa DataDomain
        data_domain = model
        model = missing
    else
        rmodel = reservoir_model(model)
        data_domain = rmodel.data_domain
    end
    cell_centroids = data_domain[:cell_centroids]
    if ismissing(aspect)
        x = cell_centroids[1, :]
        y = cell_centroids[2, :]
        xrng = maximum(x) - minimum(x)
        yrng = maximum(y) - minimum(y)
        aspect = (1.0, yrng/xrng, zaspect)
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

    if gui
        fig = plot_interactive(data_domain, arg...; z_is_depth = true, aspect = aspect, kwarg...)
        ax = fig.current_axis[]
    else
        fig, ax, plt = plot_cell_data(g, arg...; z_is_depth = true, kwarg...)
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
    for (k, w) in pairs(wells)
        tf = 0.2 + 0.1*(i/n)
        if well_color isa AbstractDict
            well_color_k = get(well_color, k, :darkred)
        else
            well_color_k = well_color
        end
        plot_well!(ax.scene, g, w;
            fontsize = well_fontsize,
            top_factor = well_top_factor_scale*tf,
            bounds_z = bounds_z,
            color = well_color_k,
            linewidth = well_linewidth,
            cell_centroids = cell_centroids,
            well_arg...
        )
        i += 1
    end
    return fig
end

# function plot_reservoir(d::DataDomain, arg...;
#         aspect = (1.0, 1.0, 1/3),
#         gui = true,
#         kwarg...
#     )
#     if gui
#         fig = plot_interactive(d, arg...; z_is_depth = true, aspect = aspect, kwarg...)
#         ax = fig.current_axis[]
#     else
#         g = physical_representation(d)
#         fig, ax, plt = plot_cell_data(g, arg...; z_is_depth = true, kwarg...)
#     end
#     return fig
# end

function plot_reservoir(case::JutulCase, arg...; kwarg...)
    if length(arg) == 0
        arg = (merge(case.parameters[:Reservoir], case.state0[:Reservoir]),)
    end
    return plot_reservoir(case.model, arg...; kwarg...)
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
