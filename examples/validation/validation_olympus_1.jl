# # The OLYMPUS benchmark model: Two-phase corner-point reservoir
# Model from the [ISAPP Optimization
# challenge](https://www.isapp2.com/optimization-challenge/reservoir-model-description.html)
#
# Two-phase complex corner-point model with primary and secondary production.
#
# For more details, see [olympus](@cite)
#
using Jutul, JutulDarcy, GLMakie, DelimitedFiles, HYPRE
using Test # hide
olympus_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("OLYMPUS_1")
case = setup_case_from_data_file(joinpath(olympus_dir, "OLYMPUS_1.DATA"), backend = :csr)
ws, states = simulate_reservoir(case, output_substates = true)
# ## Plot the reservoir
plot_reservoir(case.model, key = :porosity)
# ## Plot the saturations
plot_reservoir(case.model, states, step = 10, key = :Saturations)
# ## Load reference and set up plotting
csv_path = joinpath(olympus_dir, "REFERENCE.CSV")
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

function plot_well_comparison(response, well_names, reponse_name = "$response"; ylims = missing)
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

        val = sum(qoi.*diff([0; time_jutul]))       # hide
        val_ref = sum(qoi_ref.*diff([0; time_ref])) # hide
        @test abs(val - val_ref)/abs(val) < 0.03    # hide

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
        if !ismissing(ylims)
            ylims!(ax, ylims)
        end
    end
    l1 = LineElement(color = :black, linestyle = nothing)
    l2 = LineElement(color = :black, linestyle = :dash)

    Legend(fig[1:3, 2], linehandles, linelabels, nbanks = 3)
    Legend(fig[4, 2], [l1, l2], ["JutulDarcy.jl", "OPM Flow"])
    fig
end
# ## Plot water production rates
plot_well_comparison(:wrat, prod, "Water surface production rate")
# ## Plot oil production rates
plot_well_comparison(:orat, prod, "Oil surface production rate", ylims = (0, 5))
# ## Plot water injection rates
plot_well_comparison(:wrat, inj, "Water injection rate")
