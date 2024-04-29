

using JutulDarcy, Jutul
Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)
nx = ny = 25
nz = 10
cart_dims = (nx, ny, nz)
physical_dims = (1000.0, 1000.0, 100.0).*meter
g = CartesianMesh(cart_dims, physical_dims)
domain = reservoir_domain(g, permeability = 0.3Darcy, porosity = 0.2)
Injector = setup_vertical_well(domain, 1, 1, name = :Injector);
Producer = setup_well(domain, (nx, ny, 1), name = :Producer);


phases = (LiquidPhase(), VaporPhase())
rhoLS = 1000.0kg/meter^3
rhoGS = 100.0kg/meter^3
reference_densities = [rhoLS, rhoGS]
sys = ImmiscibleSystem(phases, reference_densities = reference_densities)

model, parameters = setup_reservoir_model(domain, sys, wells = [Injector, Producer])

"""
Machine Learning-based method for computing relative permeabilites
"""
struct MLModelRelativePermeabilities{T} <: JutulDarcy.AbstractRelativePermeabilities
    test_value::T
    function MLModelRelativePermeabilities(input_test_value)
        new{typeof(input_test_value)}(input_test_value)
    end
end

Jutul.@jutul_secondary function update_kr!(kr, kr_def::MLModelRelativePermeabilities, model, Saturations, ix)
    test_value = kr_def.test_value
    for i in ix
        for ph in axes(kr, 1)
            S = Saturations[ph, i]
            kr[ph, i] = S*test_value
        end
    end
    return kr
end

c = [1e-6, 1e-4]/bar
density = ConstantCompressibilityDensities(
    p_ref = 100*bar,
    density_ref = reference_densities,
    compressibility = c
)
#kr = BrooksCoreyRelativePermeabilities(sys, [2.0, 3.0])
#replace_variables!(model, PhaseMassDensities = density, BrooksCoreyRelativePermeabilities = kr);
kr = MLModelRelativePermeabilities(1.0)
rmodel = reservoir_model(model)
replace_variables!(rmodel, RelativePermeabilities =  kr, throw = true)

state0 = setup_reservoir_state(model,
    Pressure = 120bar,
    Saturations = [1.0, 0.0]
)
# ### Set up schedule with driving forces

nstep = 25
dt = fill(365.0day, nstep)
pv = pore_volume(model, parameters)
inj_rate = 1.5*sum(pv)/sum(dt)
rate_target = TotalRateTarget(inj_rate)
I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = rhoGS)
bhp_target = BottomHolePressureTarget(100bar)
P_ctrl = ProducerControl(bhp_target)
controls = Dict(:Injector => I_ctrl, :Producer => P_ctrl)
forces = setup_reservoir_forces(model, control = controls)

wd, states, t = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces)
##
wd(:Producer)
##
wd(:Injector, :bhp)
# ### Plot the well rates
using GLMakie

grat = wd[:Producer, :grat]
lrat = wd[:Producer, :lrat]
bhp = wd[:Injector, :bhp]

fig = Figure(size = (1200, 400))

ax = Axis(fig[1, 1],
    title = "Injector",
    xlabel = "Time / days",
    ylabel = "Bottom hole pressure / bar")
lines!(ax, t/day, bhp./bar)

ax = Axis(fig[1, 2],
    title = "Producer",
    xlabel = "Time / days",
    ylabel = "Production rate / m³/day")
lines!(ax, t/day, abs.(grat).*day)
lines!(ax, t/day, abs.(lrat).*day)

fig
# ### Launch interactive plotting of reservoir values
plot_reservoir(model, states, key = :Saturations, step = 3)
