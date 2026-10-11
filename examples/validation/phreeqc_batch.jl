# # Initial mineral inventory and equilibrium at two temperatures
# <tags: Tracers, Validation>
using Jutul, JutulDarcy, JPhreeqc, Test

mesh = CartesianMesh((2, 1), (2.0, 1.0))
model = setup_reservoir_model(reservoir_domain(mesh, porosity = 0.2),
    SinglePhaseSystem(AqueousPhase()))
rmodel = reservoir_model(model)
rmodel.parameters[:Temperature] = JutulDarcy.Temperature()
input = """
SOLUTION 1
    temp 25
    pH 9
    units mol/kgw
    Na 0.01
    Cl 0.01
    Ca 0.002
    C(4) 0.004
END
"""
setup_phreeqc_model!(model, input, equilibrium_phases = ["Calcite"], batch_size = 1)
try
    names = phreeqc_component_names(model)
    masses = phreeqc_component_molar_masses(model)
    totals = phreeqc_initial_concentrations(model)
    ca, carbon, oxygen = (findfirst(==(name), names) for name in ("Ca", "C", "O"))
    # Add 1 mmol calcite per kg aqueous carrier to conserved elemental totals.
    inventory = 1e-3
    totals[ca] += inventory*masses[ca]
    totals[carbon] += inventory*masses[carbon]
    totals[oxygen] += 3*inventory*masses[oxygen]
    state0 = setup_reservoir_state(model, Pressure = 101325.0, TracerConcentrations = totals)
    parameters = setup_parameters(rmodel, Temperature = [298.15, 323.15])
    state = Jutul.evaluate_all_secondary_variables(rmodel, state0[:Reservoir], parameters)
    aqueous = state[:PhreeqcAqueousConcentrations]
    immobile = state[:PhreeqcImmobileConcentrations]
    @test aqueous + immobile ≈ repeat(totals, 1, 2) atol = 1e-12
    @test immobile[ca, :] ./ masses[ca] ≈ immobile[carbon, :] ./ masses[carbon] rtol = 1e-6
    @test immobile[oxygen, :] ./ masses[oxygen] ≈ 3 .* immobile[carbon, :] ./ masses[carbon] rtol = 1e-6
    @test abs(aqueous[ca, 1] - aqueous[ca, 2]) > 1e-8
    minerals = state[:PhreeqcMineralConcentrations]
    @test minerals[1, :] ≈ immobile[ca, :] ./ masses[ca] rtol = 1e-12
    @test all(state[:FluidVolume] .< parameters[:StaticFluidVolume])
    occupied = state[:PhreeqcMineralVolume] .* state[:PhaseMassDensities][1, :] .* state[:FluidVolume]
    @test state[:FluidVolume] + occupied ≈ parameters[:StaticFluidVolume] rtol = 1e-12
    println("25 / 50 °C aqueous calcium (kg/kg carrier): ", aqueous[ca, :])
    println("25 / 50 °C calcite inventory (mol/kg carrier): ", immobile[ca, :] ./ masses[ca])
    println("25 / 50 °C pore volume (m³): ", state[:FluidVolume])
finally
    close_phreeqc!(model)
end
