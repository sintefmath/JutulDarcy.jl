using Jutul, JutulDarcy, HYPRE
using GLMakie
using GeoEnergyIO
data_dir = GeoEnergyIO.test_input_file_path("EGG")
data_pth = joinpath(data_dir, "EGG.DATA")
fine_case = setup_case_from_data_file(data_pth)

fine_model = fine_case.model
fine_reservoir = reservoir_domain(fine_model)
fine_mesh = physical_representation(fine_reservoir)
##
ws, states = simulate_reservoir(fine_case)
##
coarse_case = coarsen_reservoir_case(fine_case, (10, 10, 1), method = :ijk)
coarse_model = coarse_case.model
coarse_reservoir = reservoir_domain(coarse_case)
coarse_mesh = physical_representation(coarse_reservoir)

p = coarse_mesh.partition

plot_cell_data(fine_mesh, p, colormap = :devonS)
##
ws_c, states_c = simulate_reservoir(coarse_case)
##
using Statistics
wells = JutulDarcy.get_model_wells(fine_model)

p_c = states_c[end][:Pressure]
p_f = states[end][:Pressure]

caxis = extrema([extrema(p_c)..., extrema(p_f)...])

fig = Figure()
axf = Axis3(fig[1, 1], title = "Fine scale", zreversed = true)
plot_cell_data!(axf, fine_mesh, p_f, colorrange = caxis)

for (k, w) in wells
    plot_well!(axf, fine_mesh, w, fontsize = 0)
end

axc = Axis3(fig[1, 2], title = "Coarse scale", zreversed = true)
plot_cell_data!(axc, coarse_mesh, p_c, colorrange = caxis)

for (k, w) in wells
    plot_well!(axc, fine_mesh, w, fontsize = 0)
end

axf_p = Axis(fig[2, 1], ylabel = "Avg. pressure / bar")
lines!(axf_p, map(x -> mean(x[:Pressure])/1e5, states))

axf_p = Axis(fig[2, 2], ylabel = "Avg. pressure / bar")
lines!(axf_p, map(x -> mean(x[:Pressure])/1e5, states_c))

fig
##
plot_well_results([ws, ws_c], names = ["Fine", "Coarse"])
##
fine_m = reservoir_measurables(fine_case, ws, states)
coarse_m = reservoir_measurables(coarse_case, ws_c, states_c)
m = copy(fine_m)
for (k, v) in pairs(coarse_m)
    if k != :time
        m[Symbol("coarse_$k")] = v
    end
end
##
plot_reservoir_measurables(m, left = :fopr, right = :coarse_fopr, accumulated = true)