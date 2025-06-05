function setup_reservoir_model_geothermal(
        reservoir::DataDomain;
        thermal = true,
        extra_out = true,
        parameters = Dict{Symbol, Any}(),
        salt_mole_fractions = Float64[],
        salt_names = String[],
        table_arg = NamedTuple(),
        single_phase = true,
        update_reservoir = true,
        kwarg...
    )
    thermal || throw(ArgumentError("Cannot setup geothermal reservoir model with thermal = false"))
    tables_with_co2 = JutulDarcy.CO2Properties.co2_brine_property_tables(
        missing;
        salt_names = salt_names,
        salt_mole_fractions = salt_mole_fractions,
        table_arg...
    )

    function water_only_table(tabs, k)
        t = tabs[k]
        t::Jutul.BilinearInterpolant
        F = map(first, t.F)
        if k == :density
            # Density is a bit special, want to make sure that there exists some
            # derivatives in the table at low pressure and temperature since
            # otherwise the system is singular.
            ϵ = 0.0001
            for i in axes(F, 1)
                F[i, 1] *= 1.0 - ϵ
                F[i, end] *= 1.0 + ϵ
            end
            for j in axes(F, 2)
                if j == 1
                    continue
                end
                F[1, j] *= 1.0 - ϵ
                F[end, j] *= 1.0 + ϵ
            end
        end
        return Jutul.BilinearInterpolant(t.X, t.Y, F)
    end
    tables = Dict()
    for k in [:density, :heat_capacity_constant_pressure, :viscosity, :phase_conductivity]
        tables[k] = water_only_table(tables_with_co2, k)
    end

    rhoWS = first(JutulDarcy.reference_densities(:co2brine))
    if single_phase
        sys = SinglePhaseSystem(AqueousPhase(), reference_density = rhoWS)
    else
        error("Multiphase geothermal is not implemented")
    end
    # TODO: make this dynamic
    cond_water = tables[:phase_conductivity](1*si_unit(:atm), 273.15 + 20.0)
    if update_reservoir
        reservoir[:fluid_thermal_conductivity] .= cond_water
    end
    model = setup_reservoir_model(reservoir, sys; thermal = true, extra_out = false, kwarg...)
    # Tables
    rho = JutulDarcy.PressureTemperatureDependentVariable(tables[:density])
    c_p = JutulDarcy.PressureTemperatureDependentVariable(tables[:heat_capacity_constant_pressure])

    mu = JutulDarcy.PTViscosities(tables[:viscosity])

    for (k, m) in pairs(model.models)
        for (k, m) in pairs(model.models)
            if k == :Reservoir || JutulDarcy.model_or_domain_is_well(m)
                set_secondary_variables!(m;
                    PhaseMassDensities = rho,
                    PhaseViscosities = mu,
                    ComponentHeatCapacity = c_p
                )
            end
        end
    end
    rmodel = reservoir_model(model)
    outvar = rmodel.output_variables

    push!(outvar, :PhaseMassDensities)
    push!(outvar, :RockInternalEnergy)
    push!(outvar, :FluidInternalEnergy)
    push!(outvar, :TotalThermalEnergy)

    unique!(outvar)

    if extra_out
        parameters = setup_parameters(model, parameters)
        out = (model, parameters)
    else
        out = model
    end
    return out
end
