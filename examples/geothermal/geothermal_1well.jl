# # Geothermal with a single well
# This example demonstrates how to simulate a geothermal reservoir with a single
# well with a charge-discharge cycle. The reservoir has a boundary condition
# that provides a constant pressure and temperature at the boundaries.
using JutulDarcy, Jutul, HYPRE
import Dates: monthname
darcy, litre, year, second = si_units(:darcy, :litre, :year, :second)

nx = 101
nz = 100

temperature_top = convert_to_si(70.0, :Celsius)
pressure_top = convert_to_si(120.0, :bar)
temperature_surface = convert_to_si(10.0, :Celsius)

grad_p = 1000*9.81
grad_T = 0.3

# ## Set up the reservoir
g = CartesianMesh((nx, 1, nz), (250.0, 250.0, 75.0))
reservoir = reservoir_domain(g,
    permeability = [0.3, 0.3, 0.1].*darcy,
    porosity = 0.3,
    rock_thermal_conductivity = 2.0,
    fluid_thermal_conductivity = 0.6
)

depth = reservoir[:cell_centroids][3, :]
# ## Define wells and model
midx = Int(ceil(nx/2))
midz = Int(ceil(nz/2))
W = setup_vertical_well(reservoir, midx, 1, toe = midz, name = :Well)

model, parameters = setup_reservoir_model(reservoir, :geothermal, wells = W);
# ## Set up boundary and initial conditions
bcells = Int[]
pressure_res = Float64[]
temperature_res = Float64[]
for cell in 1:number_of_cells(g)
    d = depth[cell]
    push!(pressure_res, pressure_top + grad_p*d)
    push!(temperature_res, temperature_top + grad_T*d)

    I, J, K = cell_ijk(g, cell)
    if I == 1 || I == nx
        push!(bcells, cell)
    end
end

bc = flow_boundary_condition(bcells, reservoir, pressure_res[bcells], temperature_res[bcells])
# ## Set up the schedule

# ### Set up forces for charge
charge_rate = 1litre/second
discharge_rate = charge_rate

rate_target = TotalRateTarget(charge_rate)
ctrl_charge  = InjectorControl(rate_target, [1.0], density = 1000.0, temperature = temperature_surface)
forces_charge = setup_reservoir_forces(model, control = Dict(:Well => ctrl_charge), bc = bc)

# ### Set up forces for discharge
rate_target = TotalRateTarget(-discharge_rate)
ctrl_discharge = ProducerControl(rate_target)
forces_discharge = setup_reservoir_forces(model, control = Dict(:Well => ctrl_discharge), bc = bc)

# ### Set up forces for rest period
forces_rest = setup_reservoir_forces(model, bc = bc)
# ### Set up timesteps and assign forces to each timestep
num_years = 25
dt = Float64[]
forces = []
month = year/12
for year in 1:num_years
    for mno in 1:12
        mname = monthname(mno)
        if mname in ("January", "February", "March", "December")
            push!(dt, month)
            push!(forces, forces_discharge)
        elseif mname in ("June", "July", "August", "September")
            push!(dt, month)
            push!(forces, forces_charge)
        else
            @assert mname in ("April", "May", "October", "November")
            push!(dt, month)
            push!(forces, forces_rest)
        end
    end
end
# ## Set up initial state
state0 = setup_reservoir_state(model, Pressure = pressure_res, Temperature = temperature_res)
# ## Simulate the case
ws, states = simulate_reservoir(state0, model, dt,
    forces = forces,
    parameters = parameters,
    info_level = -1
)
# ## Plot the reservoir states in the interactive viewer
using GLMakie
plot_reservoir(model, states, key = :Temperature, step = num_years*12)
# ## Plot wells interactively
plot_well_results(ws)
# ## Plot the recovered energy
wd = ws.wells[:Well]
c_p_water = 4.186 # kJ/kgK

well_temp = wd[:temperature]

produced_energy = -(temperature_surface .- well_temp).*wd[:mass_rate].*c_p_water.*dt

t = cumsum(dt)./si_unit(:day)
fig = Figure()
ax = Axis(fig[1, 1], title = "Produced energy [kJ]")
lines!(ax, t, cumsum(-produced_energy))
fig
