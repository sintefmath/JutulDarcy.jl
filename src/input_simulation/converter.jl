"""
    convert_co2store_to_co2_brine(data; verbose = true)

Converts a CO2STORE data file to a co2-brine model. The conversion should be
close to equivialent for models without salt. The data will be copied if
modifications are necessary.
"""
function convert_co2store_to_co2_brine(data; verbose = true, kwarg...)
    rs = get(data, "RUNSPEC", Dict())
    if haskey(rs, "CO2STORE") || haskey(rs, "CO2STOR")
        data = convert_co2store_to_co2_brine!(deepcopy(data); verbose = verbose, kwarg...)
    else
        jutul_message("CO2STORE converter", "Model does not contain CO2STORE/CO2STOR in RUNSPEC, will not convert.", color = :green)
    end
    return data
end

"""
    convert_co2store_to_co2_brine!(data; verbose = true)

Mutating version of `convert_co2store_to_co2_brine`.
"""
function convert_co2store_to_co2_brine!(data; verbose = true)
    function cstrip(z::AbstractVector)
        z = z[1:2]
        zt = sum(z)
        @assert zt > 0 "CO2 + brine mole fractions must be above zero"
        return z./zt
    end
    function cstrip(z::AbstractMatrix)
        n = size(z, 1)
        z_new = zeros(n, 2)
        for i in 1:n
            z_new[i, :] .= cstrip(z[i, :])
        end
        return z_new
    end
    function local_msg(x)
        if verbose
            jutul_message("CO2STORE converter", x, color = :blue)
        end
    end
    function local_warn(x)
        jutul_message("CO2STORE converter", x, color = :yellow)
    end
    rs = data["RUNSPEC"]
    has_co2store = haskey(rs, "CO2STORE")
    has_co2stor = haskey(rs, "CO2STOR")
    if !(has_co2store || has_co2stor)
        throw(ArgumentError("Cannot convert if RUNSPEC is missing CO2STORE/CO2STOR keyword."))
    end
    delete!(rs, "CO2STORE")
    delete!(rs, "CO2STOR")
    delete!(rs, "COMPS")
    rs["JUTUL_CO2BRINE"] = true
    rs["OIL"] = true
    rs["GAS"] = true
    # PROPS:
    props = data["PROPS"]
    # Check CNAMES in case they are in unexpected order and warn for salts that
    # are not supported.
    cnames = get(props, "CNAMES", nothing)
    if !isnothing(cnames)
        c1 = cnames[1]
        c2 = cnames[2]
        if lowercase(c1) != "h2o"
            local_warn("First component $c1 is assumed to be H2O.")
        end
        if lowercase(c2) != "co2"
            local_warn("Second component $c2 is assumed to be CO2.")
        end
        if length(cnames) > 2
            salts = cnames[3:end]
            local_warn("More than two components were declared. Salts $(salts) will be used to adjust properties, but will not be present as individual components.")
            remap = Dict(
                "nacl" => "NaCl",
                "kcl" => "KCl",
                "caso4" => "CaSO4",
                "cacl2" => "CaCl2",
                "mgso4" => "MgSO4",
                "mgcl2" => "MgCl2"
            )
            new_salts = map(x -> remap[lowercase(x)], salts)
            if haskey(data["PROPS"], "ZMFVD")
                zmfvd = data["PROPS"]["ZMFVD"]
                if length(zmfvd) > 1
                    local_warn("Using first ZMFVD region for salts...")
                end
                if size(zmfvd, 1) > 1
                    local_warn("Using first ZMFVD entry for salts...")
                end
                zi = first(zmfvd)[1, 4:end]
                rs["SALTS"] = new_salts
                rs["SALT_MOLE_FRACTIONS"] = zi
            else
                local_warn("ZMFVD not found, will not use salts.")
            end
        end
        delete!(props, "CNAMES")
    end
    if haskey(data["PROPS"], "SALINITY")
        salinity = data["PROPS"]["SALINITY"]
        if haskey(rs, "SALTS") && salinity > 0.0
            local_warn("Both SALINITY and salts declared as COMPS are present. COMPS will be used.")
        else
            local_msg("Converting SALINITY to NaCL salt.")
            salt_mf = CO2Properties.convert_salinity_to_mole_fractions(salinity)
            data["PROPS"]["SALTS"] = ["NaCl"]
            data["PROPS"]["SALT_MOLE_FRACTIONS"] = salt_mf
        end
    end
    # ZMFVD (convert)
    zmfvd = get(props, "ZMFVD", nothing)
    if isnothing(zmfvd)
        # We can default it to be water-filled.
        props["ZMFVD"] = [[-100000.0 1.0 0.0; 100000 1.0 0.0]]
    else
        for (i, z_i) in enumerate(zmfvd)
            depth = z_i[:, 1]
            z = cstrip(z_i[:, 2:end])
            zmfvd[i] = hcat(depth, z)
        end
    end
    # Remove: MW, ACF, BIC, ZCRIT, VCRIT, TCRIT, PCRIT, EOS
    for kw in ["MW", "ACF", "BIC", "ZCRIT", "VCRIT", "TCRIT", "PCRIT", "EOS"]
        delete!(props, kw)
    end
    return data
end
