function setup_reservoir_model_co2_brine(reservoir::DataDomain,
    T = missing;
    thermal = false,
    composite = thermal,
    kwarg...
    )
    tables = co2_brine_property_tables(T)
    rho = JutulDarcy.BrineCO2MixingDensities(tables[:density])
    mu = JutulDarcy.PTViscosities(tables[:viscosity])
    if thermal
        c_v = JutulDarcy.PressureTemperatureDependentVariable(tables[:heat_capacity_constant_volume])
        rho = Pair(:flow, rho)
        mu = Pair(:flow, mu)
        c_v = Pair(:thermal, c_v)
    end
    mixture = MultiComponentMixture(["Water", "CarbonDioxide"], name = "CSP11BC-mixture")
    mixture.component_names[1] = "H₂O"
    mixture.component_names[2] = "CO₂"
    keos = KValuesEOS(tables[:K], mixture)
    # Densities
    rhoSurfaceBrine = 998.207150
    rhoSurfaceCO2 = 1.868048
    rhoS = [rhoSurfaceBrine, rhoSurfaceCO2]
    L, V = LiquidPhase(), VaporPhase()
    sys = MultiPhaseCompositionalSystemLV(keos, (L, V), reference_densities = rhoS)
    if thermal
        sys = reservoir_system(flow = sys, thermal = ThermalSystem(sys))
    elseif composite
        sys = reservoir_system(flow = sys)
    end
    model, parameters = setup_reservoir_model(reservoir, sys; reference_densities = rhoS, kwarg...);

    outvar = model[:Reservoir].output_variables
    push!(outvar, :Saturations)
    push!(outvar, :PhaseMassDensities)
    push!(outvar, :LiquidMassFractions)
    push!(outvar, :VaporMassFractions)
    if thermal
        push!(outvar, :RockInternalEnergy)
        push!(outvar, :FluidInternalEnergy)
        push!(outvar, :TotalThermalEnergy)
    end
    unique!(outvar)

    for (k, m) in pairs(model.models)
        if k == :Reservoir || JutulDarcy.model_or_domain_is_well(m)
            set_secondary_variables!(m;
                PhaseMassDensities = rho,
                Saturations = JutulDarcy.SaturationsFromDensities(),
                PhaseViscosities = mu
            )
            if thermal
                set_secondary_variables!(m;
                    ComponentHeatCapacity = c_v,
                )
            end
        end
    end
    return (model, parameters)
end

function co2_brine_property_tables(T = missing)
    getpth(n) = joinpath(@__DIR__, "tables", n)

    co2 = read_component_table(getpth("co2values.csv"))
    h2o = read_component_table(getpth("h2ovalues.csv"))
    sol = read_solubility_table(getpth("solubilities.csv"))

    prop_tab = [h2o, co2]
    K_pT = co2brine_K_values(sol, T)
    rho_pT = co2brine_phase_property_table(prop_tab, :density, T)
    mu_pT = co2brine_phase_property_table(prop_tab, :viscosity, T)
    H_pT = co2brine_phase_property_table(prop_tab, :H, T)
    C_v = co2brine_phase_property_table(prop_tab, :cv, T)
    C_p = co2brine_phase_property_table(prop_tab, :cp, T)
    E = co2brine_phase_property_table(prop_tab, :E, T)

    phase_conductivity = co2brine_phase_property_table(prop_tab, :phase_conductivity, T)

    return Dict(
        :K => K_pT,
        :density => rho_pT,
        :viscosity => mu_pT,
        :enthalpy => H_pT,
        :co2_table => co2,
        :h2o_table => h2o,
        :heat_capacity_constant_pressure => C_p,
        :heat_capacity_constant_volume => C_v,
        :internal_energy => E,
        :phase_conductivity => phase_conductivity,
        :solubility_table => sol
    )
end
