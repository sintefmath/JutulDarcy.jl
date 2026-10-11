# # PHREEQC equilibrium reactive transport in a column
# <tags: Tracers, Validation, Wells>
# Run with `julia --threads=4 --project phreeqc_column.jl` in an environment
# containing JutulDarcy and JPhreeqc. Additional species use the tracer system;
# only their aqueous split flows. Mineral amounts are encoded in total tracers.
using Jutul, JutulDarcy, JPhreeqc
using Test

nx = 12
mesh = CartesianMesh((nx, 1, 1), (12.0, 1.0, 1.0))
domain = reservoir_domain(mesh, porosity = 0.2, permeability = 1e-12)
injector = setup_well(domain, [1], name = :Injector)
producer = setup_well(domain, [nx], name = :Producer)
model = setup_reservoir_model(domain, SinglePhaseSystem(AqueousPhase()),
    wells = [injector, producer])

input = """
SOLUTION 1
    temp 25
    pH 9
    units mol/kgw
    Na 0.01
    Cl 0.01
    Ca 0.002
    C(4) 0.004
SOLUTION 0
    temp 25
    pH 7
    units mol/kgw
    Na 0.01
    Cl 0.01
    Ca 0.0001
    C(4) 0.0002
END
"""
setup_phreeqc_model!(model, input, equilibrium_phases = ["Calcite"], batch_size = 3)
rmodel = reservoir_model(model)
push!(rmodel.output_variables, :PhaseMassDensities)
try
    names = phreeqc_component_names(model)
    c0 = phreeqc_initial_concentrations(model)
    cin = phreeqc_initial_concentrations(model, solution = 0)
    ca = findfirst(==("Ca"), names)
    state0 = setup_reservoir_state(model, Pressure = 2e5, TracerConcentrations = c0)
    initial = Jutul.evaluate_all_secondary_variables(rmodel, state0[:Reservoir])
    controls = Dict(
        :Injector => InjectorControl(TotalRateTarget(1e-4), [1.0], density = 1000.0, tracers = cin),
        :Producer => ProducerControl(BottomHolePressureTarget(2e5)),
    )
    forces = setup_reservoir_forces(model, control = controls)
    result = simulate_reservoir(state0, model, fill(200.0, 20), forces = forces,
        info_level = -1, output_substates = true)
    states = result.states
    final = states[end]
    initial_mineral = initial[:PhreeqcImmobileConcentrations][ca, :]
    final_mineral = final[:PhreeqcImmobileConcentrations][ca, :]
    @test length(result.result.states) == 20
    @test all(isfinite, final[:PhreeqcAqueousConcentrations])
    @test minimum(final_mineral) >= -1e-10
    @test final_mineral[1] < initial_mineral[1]
    @test final[:FluidVolume][1] > initial[:FluidVolume][1]
    occupied = final[:PhreeqcMineralVolume] .* final[:PhaseMassDensities][1, :] .* final[:FluidVolume]
    @test final[:FluidVolume] + occupied ≈ setup_parameters(rmodel)[:StaticFluidVolume] rtol = 1e-10
    produced_ca = result.result.states[end][:Producer][:TracerConcentrations][ca, 1]
    @test abs(produced_ca - final[:PhreeqcAqueousConcentrations][ca, end]) <
        abs(produced_ca - final[:TracerConcentrations][ca, end])
    println("Initial / final immobile calcium (kg/kg carrier):")
    println(hcat(initial_mineral, final_mineral))
    println("Final aqueous calcium (kg/kg carrier):")
    println(final[:PhreeqcAqueousConcentrations][ca, :])
    println("Initial / final pore volume at the injector (m³): ",
        (initial[:FluidVolume][1], final[:FluidVolume][1]))
finally
    close_phreeqc!(model)
end
