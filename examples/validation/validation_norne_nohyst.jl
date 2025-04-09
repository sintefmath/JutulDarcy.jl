# # Norne: Real field black-oil model
# The Norne model is a real field model. The model has been adapted so that the
# input file only contains features present in JutulDarcy, with the most notable
# omissions being removal of hysteresis and threshold pressures between
# equilibriation reqgions. For more details, see the [OPM data
# webpage](https://opm-project.org/?page_id=559)
using Jutul, JutulDarcy, GLMakie, DelimitedFiles, HYPRE, GeoEnergyIO
norne_dir = GeoEnergyIO.test_input_file_path("NORNE_NOHYST")
data_pth = joinpath(norne_dir, "NORNE_NOHYST.DATA")
data = parse_data_file(data_pth)
case = setup_case_from_data_file(data);
# ## Unpack the case to see basic data structures
model = case.model
parameters = case.parameters
forces = case.forces
dt = case.dt;
# ## Plot the reservoir mesh, wells and faults
# We compose a few different plotting calls together to make a plot that shows
# the outline of the mesh, the fault structures and the well trajectories.
import Jutul: plot_mesh_edges!
reservoir = reservoir_domain(model)
mesh = physical_representation(reservoir)
wells = get_model_wells(model)
fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1], zreversed = true)
plot_mesh_edges!(ax, mesh, alpha = 0.5)
for (k, w) in wells
    plot_well!(ax, mesh, w)
end
plot_faults!(ax, mesh, alpha = 0.5)
ax.azimuth[] = -3.0
ax.elevation[] = 0.5
fig
# ## Plot the reservoir static properties in interactive viewer
fig = plot_reservoir(model, key = :porosity)
ax = fig.current_axis[]
plot_faults!(ax, mesh, alpha = 0.5)
ax.azimuth[] = -3.0
ax.elevation[] = 0.5
fig
# ## Simulate the model
ws, states = simulate_reservoir(case, output_substates = true)
# ## Plot the reservoir solution
fig = plot_reservoir(model, states, step = 247, key = :Saturations)
ax = fig.current_axis[]
ax.azimuth[] = -3.0
ax.elevation[] = 0.5
fig

# ## Load reference and set up plotting
csv_path = joinpath(norne_dir, "REFERENCE.CSV")
data_ref, header = readdlm(csv_path, ',', header = true)
time_ref = data_ref[:, 1]
time_jutul = deepcopy(ws.time)
wells = deepcopy(ws.wells)
wnames = collect(keys(wells))
nw = length(wnames)
day = si_unit(:day)
cmap = :tableau_hue_circle

inj = Symbol[]
prod = Symbol[]
for (wellname, well) in pairs(wells)
    qts = well[:wrat] + well[:orat] + well[:grat]
    if sum(qts) > 0
        push!(inj, wellname)
    else
        push!(prod, wellname)
    end
end

function plot_well_comparison(response, well_names, reponse_name = "$response"; cumulative = false)
    fig = Figure(size = (1000, 400))
    if response == :bhp
        ys = 1/si_unit(:bar)
        yl = "Bottom hole pressure / Bar"
    elseif response == :wrat
        ys = si_unit(:day)
        if cumulative
            yl = "Cumulative water rate / m³"
        else
            yl = "Water rate / m³/day"
        end
    elseif response == :grat
        ys = si_unit(:day)/1e6
        if cumulative
            yl = "Cumulative gas rate / 10⁶ m³"
        else
            yl = "Gas rate / 10⁶ m³/day"
        end
    elseif response == :orat
        ys = si_unit(:day)/(1000*si_unit(:stb))
        if cumulative
            yl = "Cumulative oil rate / 10³ stb"
        else
            yl = "Oil rate / 10³ stb/day"
        end
    else
        error("$response not ready.")
    end
    welltypes = []
    ax = Axis(fig[1:4, 1], xlabel = "Time / days", ylabel = yl)
    i = 1
    linehandles = []
    linelabels = []
    for well_name in well_names
        well = wells[well_name]
        label_in_csv = "$well_name:$response"
        ref_pos = findfirst(x -> x == label_in_csv, vec(header))
        qoi = copy(well[response]).*ys
        qoi_ref = data_ref[:, ref_pos].*ys
        tot_rate = copy(well[:rate])
        grat_ref = data_ref[:, findfirst(x -> x == "$well_name:grat", vec(header))]
        orat_ref = data_ref[:, findfirst(x -> x == "$well_name:orat", vec(header))]
        wrat_ref = data_ref[:, findfirst(x -> x == "$well_name:wrat", vec(header))]
        tot_rate_ref = grat_ref + orat_ref + wrat_ref
        if cumulative
            @. qoi_ref[tot_rate_ref == 0] = 0
            @. qoi[tot_rate == 0] = 0
            qoi_ref = cumsum(qoi_ref.*diff([0, time_ref...]./day))
            qoi = cumsum(qoi.*diff([0, time_jutul...]./day))
        else
            @. qoi_ref[tot_rate_ref == 0] = NaN
            @. qoi[tot_rate == 0] = NaN
        end
        crange = (1, max(length(well_names), 2))
        lh = lines!(ax, time_jutul./day, abs.(qoi),
            color = i,
            colorrange = crange,
            label = "$well_name", colormap = cmap
        )
        push!(linehandles, lh)
        push!(linelabels, "$well_name")
        lines!(ax, time_ref./day, abs.(qoi_ref),
            color = i,
            colorrange = crange,
            linestyle = :dash,
            colormap = cmap
        )
        i += 1
    end
    l1 = LineElement(color = :black, linestyle = nothing)
    l2 = LineElement(color = :black, linestyle = :dash)

    Legend(fig[1:3, 2], linehandles, linelabels, nbanks = 3)
    Legend(fig[4, 2], [l1, l2], ["JutulDarcy.jl", "OPM Flow"])
    fig
end
# ## Injector bhp
plot_well_comparison(:bhp, inj, "Bottom hole pressure")
# ## Gas injection rates
# ### Rates
plot_well_comparison(:grat, inj, "Gas surface injection rate")
# ## Cumulative gas injection rates
plot_well_comparison(:grat, inj, "Cumulative gas surface injection rate", cumulative = true)
# ## Water injection rates
# ### Rates
plot_well_comparison(:wrat, inj, "Water surface injection rate")
# ### Cumulative rates
plot_well_comparison(:wrat, inj, "Cumulative water surface injection rate", cumulative = true)
# ## Producer bhp
plot_well_comparison(:bhp, prod, "Bottom hole pressure")
# ## Oil production rates
# ### Rates
plot_well_comparison(:orat, prod, "Oil surface production rate")
# ### Cumulative rates
plot_well_comparison(:orat, prod, "Cumulative oil surface production rate", cumulative = true)
# ## Gas production rates
# ### Rates
plot_well_comparison(:grat, prod, "Gas surface production rate")
# ### Cumulative rates
plot_well_comparison(:grat, prod, "Cumulative gas surface production rate", cumulative = true)
# ## Water production rates
# ### Rates
plot_well_comparison(:wrat, prod, "Water surface production rate")
# ### Cumulative rates
plot_well_comparison(:wrat, prod, "Cumulative water surface production rate", cumulative = true)
# ## Interactive plotting of field statistics
plot_reservoir_measurables(case, ws, states, left = :fgpr, right = :pres)
# ## Plot wells
plot_well_results(ws)
