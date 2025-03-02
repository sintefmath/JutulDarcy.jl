# # Geothermal doublet
# This example demonstrates how to set up a geothermal doublet simulation using
# JutulDarcy. We will use two different PVT functions--one simple and one
# realistic--to highlight the importance of accurate fluid physics in geothermal
# simulations.
using Jutul, JutulDarcy, HYPRE, GeoEnergyIO, GLMakie
meter, kilogram, bar, year = si_units(:meter, :kilogram, :bar, :year)

# ## Make setup function
# We use the synthethic EGG model [egg_model](@cite) to emulate realistic
# geology. Instead of using the original wells, we set up a simple
# injector-producer doublet, placed so that injected fluids will sweep a large
# part of the reservoir. 

# Set up EGG model
egg_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("EGG")
data_pth = joinpath(egg_dir, "EGG.DATA")
case0 = setup_case_from_data_file(data_pth)
domain = reservoir_model(case0.model).data_domain;

# Make setup function
function setup_doublet(sys)

    inj_well = setup_vertical_well(domain, 45, 15, name = :Injector, simple_well = false)
    prod_well = setup_vertical_well(domain, 15, 45, name = :Producer, simple_well = false)

    model, _ = setup_reservoir_model(
        domain, sys,
        thermal = true,
        wells = [inj_well, prod_well],
    );
    rmodel = reservoir_model(model)
    push!(rmodel.output_variables, :PhaseMassDensities, :PhaseViscosities)

    state0 = setup_reservoir_state(model,
        Pressure = 50bar,
        Temperature = convert_to_si(90, :Celsius)
    )

    time = 50year
    pv_tot = sum(pore_volume(reservoir_model(model).data_domain))
    rate = 2*pv_tot/time
    rate_target = TotalRateTarget(rate)
    ctrl_inj  = InjectorControl(rate_target, [1.0],
        density = 1000.0, temperature = convert_to_si(10.0, :Celsius))

    bhp_target = BottomHolePressureTarget(25bar)
    ctrl_prod = ProducerControl(bhp_target)

    control = Dict(:Injector => ctrl_inj, :Producer => ctrl_prod)

    dt = 4year/12
    dt = fill(dt, Int(time/dt))

    forces = setup_reservoir_forces(model, control = control)

    return JutulCase(model, dt, forces, state0 = state0)

end

# ## Simple fluid physics
# We start by setting up a simple fluid physics where water is slightly
# compressible, but with no influence of temperature. Viscosity is constant.
rhoWS = 1000.0kilogram/meter^3
sys = SinglePhaseSystem(AqueousPhase(), reference_density = rhoWS)
case_simple = setup_doublet(sys)
results_simple = simulate_reservoir(case_simple);
# Interactive plot of the reservoir state
plot_reservoir(case_simple.model, results_simple.states)

# ## Realistic fluid physics
# Next, we repeat the simulation with more realistic fluid physics. We use a
# formulation from [NIST](https://webbook.nist.gov/chemistry/fluid/) where
# density, viscosity and heat capacity depend on pressure and temperature. 
case_real = setup_doublet(:geothermal)
results_real = simulate_reservoir(case_real);
# Interactive plot of the reservoir state
plot_reservoir(case_real.model, results_real.states)

# ## Compare results
# A key performace metric for geothermal doublets is the time it takes before
# the cold water injected to uphold pressure reaches the producer. At this
# point, production temperature will rapidly decline, so that the breakthrough
# time effectivelt defines the lifespan of the doublet. We plot the well results
# for the two simulations to compare the two different PVT formulations. Since
# water viscosty is not affected by temperature in the simple PVT model, water
# movement is much faster in this scenario, thereby grossly underestimating the
# lifespan of the doublet compared to the realistic PVT. This effect is further
# amplified by the thermal shrinkage due to colling present in the realistic PVT
# model.
plot_well_results([results_simple.wells, results_real.wells]; names = ["Simple", "Realistic"])

# Finally, we plot the density to see how the two simulations differ. As density
# in the the simple PVT is only dependent on pressure, it is largely constant
# except from in the vicinity of the wells, where pressure gradients are larger.
# In the realistic PVT, where density is a function of both pressure and
# temperature, we see that it is affected in all regions swept by the injected
# cold water.
ρ_simple = map(s -> s[:PhaseMassDensities], results_simple.states)
ρ_real = map(s -> s[:PhaseMassDensities], results_real.states)
Δρ = map(Δρ -> Dict(:DensityDifference => Δρ), ρ_simple .- ρ_real)
plot_reservoir(case_real.model, Δρ; step = length(Δρ))
