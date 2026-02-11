# # Group control for wells
# <tags: Immiscible, Introduction, Wells, GroupControl>

# This example demonstrates how to use group control for wells. Group control
# allows specifying a combined target for a group of wells, where each well in
# the group is allocated a fraction of the group target. This is useful for
# modelling field-level constraints where a group of producers or injectors
# share a common rate or pressure target.

# ## Preliminaries
using JutulDarcy, Jutul

Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day);

# ## Set up the reservoir domain
nx = ny = 10
nz = 4
dims = (nx, ny, nz)
g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
nlayer = nx*ny
K = vcat(
    fill(0.65, nlayer),
    fill(0.3, nlayer),
    fill(0.5, nlayer),
    fill(0.2, nlayer)
    )*Darcy

domain = reservoir_domain(g, permeability = K, porosity = 0.2)

# ## Define wells
# We create two injectors and two producers. The injectors will share a group
# rate target, and the producers will share a separate group target.
Inj1 = setup_vertical_well(domain, 1, 1, name = :Injector1)
Inj2 = setup_vertical_well(domain, 1, ny, name = :Injector2)
Prod1 = setup_vertical_well(domain, nx, 1, name = :Producer1)
Prod2 = setup_vertical_well(domain, nx, ny, name = :Producer2)

# ## Define well groups
# We define two groups: one for injectors and one for producers.
groups = Dict(
    :InjectorGroup => [:Injector1, :Injector2],
    :ProducerGroup => [:Producer1, :Producer2]
)

# ## Set up the fluid system
phases = (LiquidPhase(), VaporPhase())
rhoLS = 1000.0
rhoGS = 100.0
rhoS = [rhoLS, rhoGS] .* kg/meter^3
sys = ImmiscibleSystem(phases, reference_densities = rhoS)

# ## Create the model
# The `well_groups` keyword argument passes the group definitions to the model.
model, parameters = setup_reservoir_model(domain, sys,
    wells = [Inj1, Inj2, Prod1, Prod2],
    extra_out = true,
    well_groups = groups
)

c = [1e-6/bar, 1e-4/bar]
ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
replace_variables!(model, PhaseMassDensities = ρ);

# ## Set up initial state
state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])

# ## Set up time-steps and rates
dt = repeat([30.0]*day, 12*3)
pv = pore_volume(model, parameters)
total_inj_rate = sum(pv)/sum(dt)

# ## Define group-controlled well controls
# ### Injector group control
# The total injection rate is shared equally between both injectors. Each
# injector will target `total_inj_rate / 2` automatically via the allocation
# factor (1/nw = 1/2).
inj_group_target = TotalRateTarget(total_inj_rate)
I1_ctrl = GroupControl(
    InjectorControl(inj_group_target, [0.0, 1.0], density = rhoGS),
    :InjectorGroup
)
I2_ctrl = GroupControl(
    InjectorControl(inj_group_target, [0.0, 1.0], density = rhoGS),
    :InjectorGroup
)

# ### Producer group control
# The total production rate is shared equally between both producers. Each
# producer will target half of the group production rate.
prod_group_target = TotalRateTarget(-total_inj_rate)
P1_ctrl = GroupControl(
    ProducerControl(prod_group_target),
    :ProducerGroup
)
P2_ctrl = GroupControl(
    ProducerControl(prod_group_target),
    :ProducerGroup
)

controls = Dict(
    :Injector1 => I1_ctrl,
    :Injector2 => I2_ctrl,
    :Producer1 => P1_ctrl,
    :Producer2 => P2_ctrl
)

forces = setup_reservoir_forces(model, control = controls)

# ## Simulate the model
result = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces)
wd, states, t = result;

# ## Examine the results
# Both injectors should inject approximately equal amounts and both producers
# should produce approximately equal amounts, since the default allocation
# factor is 1/nw = 0.5 for each well in a group of 2.
x = t ./ day
println("Injector1 rate (last step): ", wd[:Injector1][:rate][end]*day, " m³/day")
println("Injector2 rate (last step): ", wd[:Injector2][:rate][end]*day, " m³/day")
println("Producer1 rate (last step): ", wd[:Producer1][:rate][end]*day, " m³/day")
println("Producer2 rate (last step): ", wd[:Producer2][:rate][end]*day, " m³/day")
println("Group injection target: ", total_inj_rate*day, " m³/day")
println("Group production target: ", -total_inj_rate*day, " m³/day")
