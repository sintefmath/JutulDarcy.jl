function read_solubility_table(name; fix = true)
    x, header = DelimitedFiles.readdlm(name, ',', header = true, skipstart = 7)
    if fix
        new_header = Symbol[]
        lookup = Dict(
            "# temperature [°C]" => :T,
            " phase pressure [Pa]" => :p,
            "         x_CO2 [-]" => :x_co2,
            "         y_H2O [-]" => :y_h2o,
        )
        for label in header
            push!(new_header, lookup[label])
        end
        header = new_header
    end
    return (data = x, header = header)
end

function read_component_table(name; fix = true)
    x, header = DelimitedFiles.readdlm(name, ',', header = true, skipstart = 7)
    if fix
        new_header = Symbol[]
        lookup = Dict(
            "# temperature [°C]" => :T,
            "     pressure [Pa]" => :p,
            "   density [kg/m3]" => :density,
            "  viscosity [Pa.s]" => :viscosity,
            "   enthalpy [J/kg]" => :H,
            " cv [J/kg.K]" => :cv,
            " cp [J/kg.K]" => :cp,
            " internal energy [J/kg.K]" => :E,
            " thermal_conductivity [W/m.K]" => :phase_conductivity
        )
        for label in header
            push!(new_header, lookup[label])
        end
        header = new_header
    end
    return (data = x, header = header)
end

function get_column(tab, x; index = false)
    ix = findfirst(isequal(x), tab.header)
    col = tab.data[:, ix]
    if index
        return (col, ix)
    else
        return col
    end
end