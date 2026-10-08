using Test, Jutul, JutulDarcy

import JutulDarcy.CO2Properties: compute_co2_brine_props

@testset "CO2-Brine Properties" begin
    p = 250e5
    T = 273.0 + 30

    props = compute_co2_brine_props(p, T, [0.0], ["NaCl"])
    @test props == compute_co2_brine_props(p, T)

    @test props[:K] ≈ [0.00447321, 37.0489] rtol = 1e-3
    @test props[:viscosity] ≈ [0.000800893, 0.000106415] rtol = 1e-3
    @test props[:density] ≈ [1006.68, 913.948] rtol = 1e-3

    p = 180e5
    T = 273.0 + 60
    props = compute_co2_brine_props(180e5, T)

    @test props[:K] ≈ [0.00843427, 46.5601] rtol = 1e-3
    @test props[:viscosity] ≈ [0.000470367, 5.96663e-5] rtol = 1e-3
    @test props[:density] ≈ [991.11, 668.644] rtol = 1e-3

    props = compute_co2_brine_props(180e5, T, [0.05], ["NaCl"])
    @test props[:K] ≈ [0.00843427, 88.0887] rtol = 1e-3
    @test props[:viscosity] ≈ [0.00064919, 5.96663e-5] rtol = 1e-3
    @test props[:density] ≈ [1092.47, 668.644] rtol = 1e-3

    props = compute_co2_brine_props(180e5, T, [0.01, 0.01, 0.005, 0.01, 0.01, 0.012], ["NaCl", "KCl", "CaSO4", "CaCl2", "MgSO4", "MgCl2"])
    @test props[:K] ≈ [0.00843427, 143.368] rtol = 1e-3
    @test props[:viscosity] ≈ [0.000501, 5.9666e-5] rtol = 1e-3
    @test props[:density] ≈ [1170.37, 668.644] rtol = 1e-3
end

@testset "CO2-Brine injector heat capacities" begin
    for (simple_well, const_p, const_T) in ((true, 101325.0, 298.15), (false, 2e7, 310.0))
        reservoir = get_1d_reservoir(2)
        well = setup_well(reservoir, [1, 2]; name = :INJ, simple_well)
        model, parameters = setup_reservoir_model(reservoir, :co2brine;
            wells = [well], thermal = true, co2_source = :csp11,
            const_p, const_T, extra_out = true)
        tables = JutulDarcy.CO2Properties.co2_brine_property_tables(co2_source = :csp11)
        expected = collect(tables[:heat_capacity_constant_volume](const_p, const_T))
        for name in (:Reservoir, :INJ)
            capacity = parameters[name][:ComponentHeatCapacity]
            @test capacity ≈ repeat(expected, 1, size(capacity, 2))
        end

        # Compare the actual well parameters with the reservoir's pure-CO2
        # energy model, rather than constructing a well state with reservoir cv.
        p, T = 3e7, 283.15
        densities = collect(tables[:density](p, T))
        state = (Pressure = [p], Saturations = reshape([0.0, 1.0], 2, 1),
            PhaseMassDensities = reshape(densities, 2, 1),
            ComponentHeatCapacity = parameters[:INJ][:ComponentHeatCapacity])
        control = InjectorControl(TotalMassRateTarget(0.035), [0.0, 1.0]; temperature = T)
        h = JutulDarcy.well_top_node_enthalpy(control, model[:INJ], state, T, 1)
        @test h ≈ parameters[:Reservoir][:ComponentHeatCapacity][2, 1]*T + p/densities[2]
    end
end
