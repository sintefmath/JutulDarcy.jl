using Test, Jutul, JutulDarcy, Adapt
using JutulDarcy.MultiComponentFlash
using JutulDarcy.StaticArrays: SVector
using JutulDarcy.ForwardDiff

function thermal_property_table(ncomp; shift = 0.0)
    pressures = [1.0e7, 2.0e7]
    temperatures = [300.0, 350.0]
    properties = [SVector{ncomp}(ntuple(c ->
        shift + 100.0*c*T + 0.001*c*p, ncomp))
        for p in pressures, T in temperatures]
    return Jutul.BilinearInterpolant(pressures, temperatures, properties)
end

@testset "Geothermal steam-table import" begin
    domain = get_1d_reservoir(2)
    system = SinglePhaseSystem(AqueousPhase())
    fluid_model = Dict(
        "WaterCompressibilities" => Dict("RefPressure" => 1.0e5, "Viscosity" => 1.0e-3),
        "WaterThermalModel" => Dict(
            "LiquidDensityModel" => "STEAM_TABLE",
            "LiquidEnthalpyModel" => "STEAM_TABLE",
            "LiquidViscosityModel" => "STEAM_TABLE"))
    variables = JutulDarcy.AFISetup.setup_pvt_variables_single_phase_water(
        nothing, system, domain, fluid_model)
    @test variables[:FluidEnthalpy] isa JutulDarcy.PressureTemperatureDependentEnthalpy
    @test variables[:FluidInternalEnergy] isa JutulDarcy.FluidInternalEnergyFromFluidEnthalpy
    @test !haskey(variables, :ComponentHeatCapacity)
end

function thermal_property_model(system)
    domain = get_1d_reservoir(2)
    return SimulationModel(discretized_domain_tpfv_flow(domain), system;
        data_domain = domain)
end

@testset "Tabulated phase enthalpy and internal energy" begin
    table = thermal_property_table(2)
    variable = JutulDarcy.PressureTemperatureDependentEnthalpy(
        (table, thermal_property_table(2; shift = 1000.0)); regions = [1, 2])
    energy = JutulDarcy.FluidInternalEnergyFromFluidEnthalpy()
    system = ImmiscibleSystem((LiquidPhase(), VaporPhase()))
    model = thermal_property_model(system)
    state = (Pressure = [1.2e7, 1.8e7], Temperature = [315.0, 335.0],
        PhaseMassDensities = [1000.0 990.0; 500.0 600.0])
    H = zeros(2, 2)
    Jutul.update_secondary_variable!(H, variable, model, state, 1:2)
    for c in 1:2
        @test H[:, c] ≈ variable.tab[c](state.Pressure[c], state.Temperature[c])
    end
    U = similar(H)
    Jutul.update_secondary_variable!(U, energy, model,
        merge(state, (; FluidEnthalpy = H)), 1:2)
    @test U ≈ H .- state.Pressure' ./ state.PhaseMassDensities
    @test Jutul.get_dependencies(variable, model) == (:Pressure, :Temperature)
    @test Jutul.get_dependencies(energy, model) ==
        (:FluidEnthalpy, :Pressure, :PhaseMassDensities)

    fill!(H, -1.0)
    Jutul.update_secondary_variable!(H, variable, model, state, 2:2)
    @test H[:, 1] == [-1.0, -1.0]
    @test H[:, 2] ≈ variable.tab[2](state.Pressure[2], state.Temperature[2])

    subvariable = Jutul.subvariable(variable, Jutul.FiniteVolumeGlobalMap([2, 1], Int[]))
    @test subvariable.regions == [2, 1]
    @test subvariable.tab === variable.tab

    # A scalar table also represents a single phase, as used by geothermal.
    single_phase = thermal_property_model(SinglePhaseSystem(AqueousPhase()))
    scalar_table = Jutul.BilinearInterpolant(table.X, table.Y, first.(table.F))
    scalar_variable = JutulDarcy.PressureTemperatureDependentEnthalpy(scalar_table)
    H_single = zeros(1, 2)
    Jutul.update_secondary_variable!(H_single, scalar_variable, single_phase, state, 1:2)
    @test vec(H_single) ≈ [scalar_table(state.Pressure[c], state.Temperature[c]) for c in 1:2]
end

@testset "Compositional enthalpy mixing and derivatives" begin
    mixture = MultiComponentMixture(["Water", "CarbonDioxide", "Methane"])
    eos = KValuesEOS((p, T) -> SVector(0.01, 10.0, 20.0), mixture)
    table = thermal_property_table(3)
    variable = JutulDarcy.PressureTemperatureDependentEnthalpy(table)
    system = MultiPhaseCompositionalSystemLV(eos, (VaporPhase(), LiquidPhase()))
    model = thermal_property_model(system)
    l, v = JutulDarcy.phase_indices(system)
    state = (Pressure = [1.2e7, 1.8e7], Temperature = [315.0, 335.0],
        LiquidMassFractions = [0.8 0.6; 0.15 0.3; 0.05 0.1],
        VaporMassFractions = [0.1 0.2; 0.6 0.3; 0.3 0.5])
    @test Jutul.values_per_entity(model, variable) == 2
    @test Jutul.get_dependencies(variable, model) ==
        (:Pressure, :Temperature, :LiquidMassFractions, :VaporMassFractions)
    H = zeros(2, 2)
    Jutul.update_secondary_variable!(H, variable, model, state, 1:2)
    for c in 1:2
        component_H = table(state.Pressure[c], state.Temperature[c])
        @test H[l, c] ≈ sum(component_H .* state.LiquidMassFractions[:, c])
        @test H[v, c] ≈ sum(component_H .* state.VaporMassFractions[:, c])
    end

    energy = JutulDarcy.FluidInternalEnergyFromFluidEnthalpy()
    function liquid_energy(pT)
        p, T = pT
        rho = 1000.0 + 1e-6*p - 0.2*(T - 300.0)
        local_state = merge(state, (Pressure = [p], Temperature = [T],
            PhaseMassDensities = fill(rho, 2, 1)))
        h = similar(pT, 2, 1)
        u = similar(h)
        Jutul.update_secondary_variable!(h, variable, model, local_state, 1:1)
        Jutul.update_secondary_variable!(u, energy, model,
            merge(local_state, (; FluidEnthalpy = h)), 1:1)
        return u[l, 1]
    end
    p, T = state.Pressure[1], state.Temperature[1]
    rho = 1000.0 + 1e-6*p - 0.2*(T - 300.0)
    weighted_component = sum((1:3) .* state.LiquidMassFractions[:, 1])
    expected = [0.001*weighted_component - 1/rho + p/rho^2*1e-6,
        100.0*weighted_component - p/rho^2*0.2]
    @test ForwardDiff.gradient(liquid_energy, [p, T]) ≈ expected

    context = KernelAbstractionsContext(Jutul.KernelExecution.KernelAbstractions.CPU();
        float_type = Float32, index_type = Int32)
    adapted = Adapt.adapt(context, variable)
    @test eltype(adapted.tab[1].F) === SVector{3, Float32}
    @test Adapt.adapt(context, energy) === energy
    float_state = map(x -> Float32.(x), state)
    float_H = Matrix{Any}(undef, 2, 2)
    Jutul.update_secondary_variable!(float_H, adapted, model, float_state, 1:2)
    @test all(h isa Float32 for h in float_H)
    @test float_H ≈ H rtol = 1e-6
end
