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

@testset "CO2-Brine tabulated thermal properties" begin
    tables = JutulDarcy.CO2Properties.co2_brine_property_tables(co2_source = :csp11)
    for simple_well in (true, false), co2_physics in (:kvalue, :immiscible)
        reservoir = get_1d_reservoir(2)
        well = setup_well(reservoir, [1, 2]; name = :INJ, simple_well)
        model, parameters = setup_reservoir_model(reservoir, :co2brine;
            wells = [well], thermal = true, co2_source = :csp11,
            co2_physics, extra_out = true)
        composition = co2_physics == :kvalue ?
            (OverallMoleFractions = [0.8, 0.2],) : (Saturations = [0.8, 0.2],)
        state0 = setup_reservoir_state(model; Pressure = 3e7, Temperature = 320.0,
            composition...)
        for name in (:Reservoir, :INJ)
            m = model[name]
            variables = m.secondary_variables
            @test variables[:FluidEnthalpy] isa JutulDarcy.PressureTemperatureDependentEnthalpy
            @test variables[:FluidInternalEnergy] isa JutulDarcy.FluidInternalEnergyFromFluidEnthalpy
            @test !haskey(parameters[name], :ComponentHeatCapacity)
            state = Jutul.evaluate_all_secondary_variables(m, state0[name], parameters[name])
            for c in eachindex(state[:Pressure])
                pure_H = tables[:enthalpy](state[:Pressure][c], state[:Temperature][c])
                if co2_physics == :kvalue
                    l, v = JutulDarcy.phase_indices(m.system)
                    @test state[:FluidEnthalpy][l, c] ≈ sum(pure_H .* state[:LiquidMassFractions][:, c])
                    @test state[:FluidEnthalpy][v, c] ≈ sum(pure_H .* state[:VaporMassFractions][:, c])
                else
                    @test state[:FluidEnthalpy][:, c] ≈ pure_H
                end
                @test state[:FluidInternalEnergy][:, c] ≈
                    state[:FluidEnthalpy][:, c] .- state[:Pressure][c] ./ state[:PhaseMassDensities][:, c]
            end
        end

        # Injector temperature and mixture can differ from the local well state.
        p, T = 3e7, 283.15
        state = (Pressure = [p],)
        evaluator = JutulDarcy.injection_enthalpy_evaluator(model[:INJ])
        for mix in ([1.0, 0.0], [0.0, 1.0], [0.25, 0.75])
            control = InjectorControl(TotalMassRateTarget(0.035), mix; temperature = T)
            h = JutulDarcy.well_top_node_enthalpy(control, model[:INJ], state, T, 1, evaluator)
            @test h ≈ sum(mix .* tables[:enthalpy](p, T))
        end
        for enthalpy in (12345.0, (p, T) -> p/T)
            control = InjectorControl(TotalMassRateTarget(0.035), [0.0, 1.0]; temperature = T, enthalpy)
            h = JutulDarcy.well_top_node_enthalpy(control, model[:INJ], state, T, 1, evaluator)
            @test h ≈ (enthalpy isa Real ? enthalpy : enthalpy(p, T))
        end
    end
end

@testset "Geothermal tabulated thermal properties" begin
    tables = JutulDarcy.Geothermal.geothermal_setup_tables()
    for simple_well in (true, false)
        reservoir = get_1d_reservoir(2)
        well = setup_well(reservoir, [1, 2]; name = :INJ, simple_well)
        model, parameters = setup_reservoir_model(reservoir, :geothermal;
            wells = [well], extra_out = true)
        state0 = setup_reservoir_state(model; Pressure = 2e7, Temperature = 330.0)
        for name in (:Reservoir, :INJ)
            m = model[name]
            @test m.secondary_variables[:FluidEnthalpy] isa JutulDarcy.PressureTemperatureDependentEnthalpy
            @test m.secondary_variables[:FluidInternalEnergy] isa JutulDarcy.FluidInternalEnergyFromFluidEnthalpy
            state = Jutul.evaluate_all_secondary_variables(m, state0[name], parameters[name])
            for c in eachindex(state[:Pressure])
                p, T = state[:Pressure][c], state[:Temperature][c]
                h = tables[:enthalpy](p, T)
                @test state[:FluidEnthalpy][1, c] ≈ h
                @test state[:FluidInternalEnergy][1, c] ≈ h - p/state[:PhaseMassDensities][1, c]
            end
        end
        p, T = 2e7, 310.0
        control = InjectorControl(TotalMassRateTarget(0.035), [1.0]; temperature = T)
        evaluator = JutulDarcy.injection_enthalpy_evaluator(model[:INJ])
        h = JutulDarcy.well_top_node_enthalpy(control, model[:INJ], (; Pressure = [p]), T, 1, evaluator)
        @test h ≈ tables[:enthalpy](p, T)
    end
end
