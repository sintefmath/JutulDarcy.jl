
function plot_well!

end

function plot_well_results

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

function plot_reservoir(model, arg...; well_fontsize = 18, well_linewidth = 3, kwarg...)
    rmodel = reservoir_model(model)
    data_domain = rmodel.data_domain
    cell_centroids = data_domain[:cell_centroids]
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
    fig = plot_interactive(data_domain, arg...; kwarg...)
    g = physical_representation(data_domain)
    ax = fig.current_axis[]
    wells = Dict{Symbol, Any}()
    for (k, m) in pairs(model.models)
        w = physical_representation(m.data_domain)
        if w isa WellDomain
            wells[k] = w
        end
    end

    i = 1
    n = length(wells)
    for (k, w) in pairs(wells)
        tf =  0.2 + 0.1*(i/n)
        plot_well!(ax.scene, g, w,
            fontsize = well_fontsize,
            top_factor = tf,
            bounds_z = bounds_z,
            linewidth = well_linewidth,
            cell_centroids = cell_centroids)
        i += 1
    end
    return fig
end

function simulate_reservoir_parray(case, mode = :mpi; kwarg...)
    sim, cfg = setup_reservoir_simulator(case; mode = mode, kwarg...)
    return simulate!(sim, case.dt, forces = case.forces, config = cfg)
end

function setup_reservoir_simulator_parray

end
