
# ## Preliminarie
using JutulDarcy, Jutul

Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day);

nx = ny = 5
nz = 4
dims = (nx, ny, nz)
g = reservoir_mesh(dims, (2000.0, 1500.0, 50.0))

domain = reservoir_domain(g, permeability = 0.65*Darcy, porosity = 0.2)

Prod = setup_vertical_well(domain, 1, 1, name = :P);

I1 = setup_well(domain, [(nx, ny, 1)], name = :I1);
I2 = setup_well(domain, [(1, ny, 1)], name = :I2);
phases = (LiquidPhase(), VaporPhase())
rhoLS = 1000.0
rhoGS = 100.0
rhoS = [rhoLS, rhoGS] .* kg/meter^3
sys = ImmiscibleSystem(phases, reference_densities = rhoS)

groups = Dict(
    :G1 => [:I1, :I2]
)
model, parameters = setup_reservoir_model(domain, sys,
    wells = [I1, I2, Prod],
    well_groups = groups,
    extra_out = true
)

# This concludes the setup of the model.
state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])

dt = repeat([30.0]*day, 12*5)
pv = pore_volume(model, parameters)
inj_rate = sum(pv)/sum(dt)

I_ctrl = setup_injector_control(inj_rate/2, :rate, [0.0, 1.0], density = rhoGS)
I_ctrl = setup_injector_control(inj_rate, :group_rate, [0.0, 1.0], density = rhoGS)
P_ctrl = setup_producer_control(50*bar, :bhp)

# rate_target = TotalRateTarget(inj_rate)
# I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = rhoGS)
# bhp_target = BottomHolePressureTarget(50*bar)
# P_ctrl = ProducerControl(bhp_target)
ENV["JULIA_DEBUG"] = "JutulDarcy"
controls = Dict()
controls[:I1] = I_ctrl
controls[:I2] = I_ctrl
controls[:P] = P_ctrl
forces = setup_reservoir_forces(model, control = controls)
result = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces);
