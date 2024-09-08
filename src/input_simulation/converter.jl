"""
    convert_co2store_to_co2_brine(data; verbose = true)

Converts a CO2STORE data file to a co2-brine model. The conversion should be
close to equivialent for models without salt. The data will be copied if
modifications are necessary.
"""
function convert_co2store_to_co2_brine(data; verbose = true, kwarg...)
    if haskey(data, "RUNSPEC") && haskey(data["RUNSPEC"], "CO2STORE")
        data = convert_co2store_to_co2_brine!(deepcopy(data); verbose = verbose, kwarg...)
    else
        jutul_message("CO2STORE converter", "Model does not contain CO2STORE in RUNSPEC, will not convert.", color = :green)
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
    haskey(rs, "CO2STORE") || throw(ArgumentError("Cannot convert if RUNSPEC is missing CO2STORE keyword."))
    delete!(rs, "CO2STORE")
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
            local_warn("More than two components were declared. Salts $(cnames[3:end]) will be ignored.")
        end
        delete!(props, "CNAMES")
    end
    # ZMFVD (convert)
    zmfvd = get(props, "ZMFVD", nothing)
    if !isnothing(zmfvd)
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
