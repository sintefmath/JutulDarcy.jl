function setup_reservoir_model_co2_brine(reservoir::DataDomain;
        temperature = missing,
        thermal = false,
        composite = thermal,
        extra_out = true,
        kwarg...
    )
    tables = co2_brine_property_tables(temperature)
    rho = JutulDarcy.BrineCO2MixingDensities(tables[:density])
    mu = JutulDarcy.PTViscosities(tables[:viscosity])
    if thermal
        c_v = JutulDarcy.PressureTemperatureDependentVariable(tables[:heat_capacity_constant_volume])
    end
    mixture = MultiComponentMixture(["Water", "CarbonDioxide"], name = "CSP11BC-mixture")
    mixture.component_names[1] = "H2O"
    mixture.component_names[2] = "CO2"
    keos = KValuesEOS(tables[:K], mixture)
    # Densities
    rhoS = JutulDarcy.reference_densities(:co2brine)
    phases = JutulDarcy.get_phases(:co2brine)
    sys = MultiPhaseCompositionalSystemLV(keos, phases, reference_densities = rhoS)

    tmp = setup_reservoir_model(reservoir, sys; thermal = thermal, kwarg...);
    out = setup_reservoir_model(reservoir, sys; thermal = thermal, extra_out = extra_out, kwarg...)
    if extra_out
        model = out[1]
    else
        model = out
    end

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
    return out
end

function JutulDarcy.reference_densities(::Val{:co2brine})
    rhoSurfaceBrine = 998.207150
    rhoSurfaceCO2 = 1.868048
    return (rhoSurfaceBrine, rhoSurfaceCO2)
end

function JutulDarcy.get_phases(::Val{:co2brine})
    return (LiquidPhase(), VaporPhase())
end


function co2_brine_property_tables(T = missing; basepath = joinpath(artifact"CO2Tables_CSP11", "csp11"))
    ispath(basepath) || throw(ArgumentError("basepath $basepath does not exist."))
    getpth(n) = joinpath(basepath, n)

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

function JutulDarcy.select_injector_mixture_spec(sys::Val{:co2brine}, name, streams, type)
    rho_w, rho_co2 = JutulDarcy.reference_densities(sys)
    if lowercase(type) == "gas"
        rho = rho_co2
        mix = [0.0, 1.0]
        phases_mix = ((1, 0.0), (2, 1.0))
    else
        rho = rho_w
        mix = [1.0, 0.0]
        phases_mix = ((1, 1.0), (2, 0.0))
    end
    return (rho, mix, phases_mix)
end
