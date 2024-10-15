function setup_reservoir_model_co2_brine(reservoir::DataDomain;
        temperature = missing,
        thermal = false,
        co2_physics = :kvalue,
        co2_table_directory = missing,
        co2_source = missing,
        co2_density = :nist,
        extra_out = true,
        salt_names = String[],
        salt_mole_fractions = Float64[],
        kwarg...
    )
    tables = co2_brine_property_tables(temperature,
        basepath = co2_table_directory,
        salt_names = salt_names,
        salt_mole_fractions = salt_mole_fractions,
        co2_source = co2_source,
        co2_density = co2_density
    )
    rho = JutulDarcy.BrineCO2MixingDensities(tables[:density])
    mu = JutulDarcy.PTViscosities(tables[:viscosity])
    if thermal
        c_v = JutulDarcy.PressureTemperatureDependentVariable(tables[:heat_capacity_constant_volume])
    end
    rhoS = JutulDarcy.reference_densities(:co2brine)
    phases = JutulDarcy.get_phases(:co2brine)

    if co2_physics == :kvalue
        is_compositional = true
        mixture = MultiComponentMixture(["Water", "CarbonDioxide"], name = "CSP11BC-mixture")
        mixture.component_names[1] = "H2O"
        mixture.component_names[2] = "CO2"
        keos = KValuesEOS(tables[:K], mixture)
        # Densities
        sys = MultiPhaseCompositionalSystemLV(keos, phases, reference_densities = rhoS)
    elseif co2_physics == :immiscible
        is_compositional = false
        sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    else
        error("Unknown physics argument for co2_physics: $co2_physics is not one of :kvalue or :immiscible.")
    end

    out = setup_reservoir_model(reservoir, sys; thermal = thermal, extra_out = extra_out, kwarg...)
    if extra_out
        model = out[1]
    else
        model = out
    end

    outvar = model[:Reservoir].output_variables
    push!(outvar, :Saturations)
    push!(outvar, :PhaseMassDensities)
    if is_compositional
        push!(outvar, :LiquidMassFractions)
        push!(outvar, :VaporMassFractions)
    end
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
                PhaseViscosities = mu
            )
            if is_compositional
                set_secondary_variables!(m;
                    Saturations = JutulDarcy.SaturationsFromDensities()
                )
            end
            if thermal
                set_secondary_variables!(m;
                    ComponentHeatCapacity = c_v,
                )
            elseif !is_compositional
                set_parameters!(m, Temperature = JutulDarcy.Temperature())
            end
        end
    end
    return out
end

function replace_co2_brine_properties!(dens, visc, K, salt_mole_fractions, salt_names; co2_density::Symbol, T_value)
    p = dens.X
    if ismissing(T_value)
        T = dens.Y
    else
        T = [T_value]
    end
    co2_density in (:rk, :nist) || throw(ArgumentError("co2_density argument must be :rk (for Redlich-Kwong) or :nist (for NIST tables)"))
    replace_co2 = co2_density != :nist
    for (i, p_i) in enumerate(p)
        if i == 1
            p_i = p[2]
        elseif i == length(p)
            p_i = p[end-1]
        end
        for (j, T_i) in enumerate(T)
            if length(T) == 1
                T_i = only(T)
            else
                if j == 1
                    T_i = T[2]
                elseif j == length(T)
                    T_i = T[end-1]
                end
            end
            props = compute_co2_brine_props(p_i, T_i, salt_mole_fractions, salt_names, check = false)
            dens_calc = props[:density]
            if replace_co2
                dens.F[i, j] = dens_calc
            else
                pair_type = eltype(dens.F)
                v_co2 = dens.F[i, j][2]
                dens.F[i, j] = pair_type(dens_calc[1], v_co2)
            end
            dens.F[i, j] = props[:density]
            visc.F[i, j] = props[:viscosity]
            K.F[i, j] = props[:K]
        end
    end
end

function JutulDarcy.reference_densities(::Val{:co2brine})
    rhoSurfaceBrine = 998.207150
    rhoSurfaceCO2 = 1.868048
    return (rhoSurfaceBrine, rhoSurfaceCO2)
end

function JutulDarcy.get_phases(::Val{:co2brine})
    return (LiquidPhase(), VaporPhase())
end

"""
    co2_brine_property_tables(T = missing;
        basepath = missing,
        salt_names = salt_names,
        salt_mole_fractions = salt_mole_fractions,
        co2_source = missing,
        co2_density = :rk
    )

Set up brine-CO2 property tables for a wide range of p, T. If a single
temperature value `T` is provided, an isothermal table will be generated for
that temperature.


Keyword arguments:

    - `basepath`: Path to folder containing tables. This can be used to override the default tables.
    - `salt_names`: Names of salts to in include (passed onto `compute_co2_brine_props`)
    - `salt_mole_fractions`: Corresponding salt mole fractions (same length as `salt_names`)
    - `co2_source = missing`: Source of CO2 values. Can be:
        - `:table`: Use table directly.
        - `:salo24`: Use correlations implemented in paper by Salo et al (support for salts)
        - `:csp11`: Use values exactly as given in CSP11 benchmark (no salts)
    - `co2_density=:rk`: How to obtain density of pure CO2 phase. Can be `:nist` for NIST tables or `:rk` for adjusted Redlich-Kwong from Spycher paper.

"""
function co2_brine_property_tables(T = missing;
        basepath = missing,
        salt_names = String[],
        salt_mole_fractions = Float64[],
        co2_source = missing,
        co2_density = :rk
    )
    if ismissing(co2_source)
        co2_source = ifelse(ismissing(basepath), :salo24, :table)
    end
    if ismissing(basepath)
        basepath = joinpath(artifact"CO2Tables_CSP11", "csp11")
    end
    length(salt_mole_fractions) == length(salt_names) || throw(ArgumentError("salt_names ($salt_names) and salt_mole_fractions ($salt_mole_fractions) must have equal length."))

    ispath(basepath) || throw(ArgumentError("basepath $basepath does not exist."))
    getpth(n) = joinpath(basepath, n)

    co2 = read_component_table(getpth("co2values.csv"))
    h2o = read_component_table(getpth("h2ovalues.csv"))
    sol = read_solubility_table(getpth("solubilities.csv"))

    prop_tab = [h2o, co2]
    K_pT = co2brine_K_values(sol, T)
    get_tab(k::Symbol) = co2brine_phase_property_table(prop_tab, k, T)
    rho_pT = get_tab(:density)
    mu_pT = get_tab(:viscosity)
    H_pT = get_tab(:H)
    C_v = get_tab(:cv)
    C_p = get_tab(:cp)
    E = get_tab(:E)
    phase_conductivity = get_tab(:phase_conductivity)

    if co2_source == :table
        if length(salt_mole_fractions) > 0
            jutul_message("CO2-Brine", "Salts $salt_names were provided but table was also provided as $basepath (or use of table without modificatiosn was specified by setting co2_source = :table). Table will not be adjusted for salinity.", color = :yellow)
        end
    elseif co2_source == :salo24 || length(salt_names) > 0
        replace_co2_brine_properties!(rho_pT, mu_pT, K_pT.K, salt_mole_fractions, salt_names, co2_density = co2_density, T_value = T)
    else
        co2_source == :csp11 || throw(ArgumentError("co2_source argument must be either :csp11, :table or :salo24"))
    end

    tables = Dict(
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
    return tables
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
