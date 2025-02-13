# # The Egg model: Two-phase oil-water model
# A two-phase model that is taken from the first member of the EGG ensemble. The
# model is a synthetic case with channelized permeability and water injection
# with fixed controls. For more details, see the paper where the ensemble is
# introduced:
#
# [Jansen, Jan-Dirk, et al. "The egg model–a geological ensemble for reservoir
# simulation." Geoscience Data Journal 1.2 (2014): 192-195.](https://doi.org/10.1002/gdj3.21)
using Jutul, JutulDarcy, GLMakie, DelimitedFiles, HYPRE
egg_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("EGG")
case = setup_case_from_data_file(joinpath(egg_dir, "EGG.DATA"))
ws, states = simulate_reservoir(case, output_substates = true)
# ## Plot the reservoir solution
plot_reservoir(case.model, states, step = 135, key = :Saturations)

# ## Load reference solution (OPM Flow)
# We load a CSV file with the reference solution and set up plotting
csv_path = joinpath(egg_dir, "REFERENCE.CSV")
data, header = readdlm(csv_path, ',', header = true)
time_ref = data[:, 1]
time_jutul = deepcopy(ws.time)
wells = deepcopy(ws.wells)
wnames = collect(keys(wells))
nw = length(wnames)
day = si_unit(:day)
cmap = :tableau_hue_circle

inj = Symbol[]
prod = Symbol[]
for (wellname, well) in pairs(wells)
    qts = well[:wrat] + well[:orat]
    if sum(qts) > 0
        push!(inj, wellname)
    else
        push!(prod, wellname)
    end
end

function plot_well_comparison(response, well_names, reponse_name = "$response")
    fig = Figure(size = (1000, 400))
    if response == :bhp
        ys = 1/si_unit(:bar)
        yl = "Bottom hole pressure / Bar"
    elseif response == :wrat
        ys = si_unit(:day)
        yl = "Surface water rate / m³/day"
    elseif response == :orat
        ys = si_unit(:day)/(1000*si_unit(:stb))
        yl = "Surface oil rate / 10³ stb/day"
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
        qoi_ref = data[:, ref_pos].*ys
        tot_rate = copy(well[:rate])
        @. qoi[tot_rate == 0] = NaN
        orat_ref = data[:, findfirst(x -> x == "$well_name:orat", vec(header))]
        wrat_ref = data[:, findfirst(x -> x == "$well_name:wrat", vec(header))]
        tot_rate_ref = orat_ref + wrat_ref
        @. qoi_ref[tot_rate_ref == 0] = NaN
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
    Legend(fig[4, 2], [l1, l2], ["JutulDarcy.jl", "E100"])
    fig
end
# ## Well responses and comparison
# As the case is a two-phase model with water injection, we limit the results to
# plots of the producer water and oil rates.

# ### Water production rates
plot_well_comparison(:wrat, prod, "Producer water surface rate")

# ### Oil production rates
plot_well_comparison(:orat, prod, "Producer oil surface rate")
