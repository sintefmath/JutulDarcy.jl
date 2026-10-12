# # Water injection with the same fracture objects in EFM and cut-cell DFM
using Jutul, JutulDarcy

matrix = reservoir_domain(reservoir_mesh((8, 8, 4), (80.0, 80.0, 20.0));
    permeability = 0.1si_unit(:darcy), porosity = 0.2)
fractures = [
    RectangleFracture([40, 40, 10], [1, 0, 0], (80, 20); aperture = 1e-3,
        permeability = 100si_unit(:darcy)),
    EllipticalFracture([40, 40, 10], [0, 1, 0], (30, 8); direction = [1, 0, 0],
        segments = 16, aperture = 1e-3, permeability = 100si_unit(:darcy))]
wells = [setup_vertical_well(matrix, 1, 1; name = :Injector, simple_well = true),
    setup_vertical_well(matrix, 8, 8; name = :Producer, simple_well = true)]
system = ImmiscibleSystem((AqueousPhase(), LiquidPhase());
    reference_densities = (1000.0, 1000.0))

efm = setup_fractured_reservoir_model(matrix, fractures, system; method = :efm, wells = wells)
dfm = setup_fractured_reservoir_model(matrix, fractures, system; method = :dfm, wells = wells)

control = Dict(:Injector => InjectorControl(TotalRateTarget(10/si_unit(:day)), [1.0, 0.0]; density = 1000.0),
    :Producer => ProducerControl(BottomHolePressureTarget(100si_unit(:bar))))
dt = fill(10si_unit(:day), 12)
cases = Dict{Symbol, JutulCase}()
results = Dict{Symbol, Any}()
for (name, model) in ((:EFM, efm), (:DFM, dfm))
    replace_variables!(model, RelativePermeabilities = BrooksCoreyRelativePermeabilities(system, 1.0))
    parameters = setup_parameters(model)
    for submodel in (:Reservoir, :Fractures, :Injector, :Producer)
        parameters[submodel][:PhaseViscosities] .= 1e-3
    end
    state0 = setup_reservoir_state(model; Pressure = 100si_unit(:bar), Saturations = [0.0, 1.0])
    forces = setup_reservoir_forces(model; control = control)
    case = JutulCase(model, dt, forces; state0 = state0, parameters = parameters)
    cases[name] = case
    results[name] = simulate_reservoir(case)
end

# Load GLMakie to inspect either result:
# plot_reservoir(cases[:EFM], results[:EFM]; key = :Saturations)
# plot_reservoir(cases[:DFM], results[:DFM]; key = :Saturations)
