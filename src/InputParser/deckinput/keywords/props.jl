
function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:RPTPROPS})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FULLIMP})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:CNAMES})
    n = compositional_number_of_components(outer_data)
    templ = fill("", n)
    rec = read_record(f)
    data["CNAMES"] = parse_defaulted_line(rec, templ)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:EOS})
    rec = read_record(f)
    data["EOS"] = only(rec)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:STCOND})
    std = parse_deck_vector(f)
    @assert length(std) == 2
    swap_unit_system_axes!(std, units, [:relative_temperature, :pressure])
    data["STCOND"] = std
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:BIC})
    n = compositional_number_of_components(outer_data)
    bic = parse_deck_vector(f)
    @assert length(bic) == n*(n-1)÷2 "Bad length for BIC input."
    m = zeros(n, n)
    ix = 1
    for i in 1:n
        for j in 1:(i-1)
            m[i, j] = bic[ix]
            ix += 1
        end
    end
    data["BIC"] = Symmetric(collect(m'))
end

function parse_compositional_helper!(f, outer_data, data, k)
    n = compositional_number_of_components(outer_data)
    val = parse_deck_vector(f)
    @assert length(val) == n "One $k should be provided per component (expected $n, was $(length(val)))."
    return val
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ACF})
    data["ACF"] = parse_compositional_helper!(f, outer_data, data, "ACF")
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PCRIT})
    p_c = parse_compositional_helper!(f, outer_data, data, "PCRIT")
    swap_unit_system!(p_c, units, :pressure)
    data["PCRIT"] = p_c
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:TCRIT})
    p_c = parse_compositional_helper!(f, outer_data, data, "TCRIT")
    swap_unit_system!(p_c, units, :absolute_temperature)
    data["TCRIT"] = p_c
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MW})
    mw = parse_compositional_helper!(f, outer_data, data, "MW")
    swap_unit_system!(mw, units, :molar_mass)
    data["MW"] = mw
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:VCRIT})
    V = parse_compositional_helper!(f, outer_data, data, "VCRIT")
    swap_unit_system!(V, units, :critical_volume)
    data["VCRIT"] = V
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:SWOF}, Val{:SGOF}})
    k = unpack_val(v)
    sat_tab = parse_saturation_table(f, outer_data)
    for tab in sat_tab
        swap_unit_system_axes!(tab, units, (:identity, :identity, :identity, :pressure))
    end
    data["$k"] = sat_tab
end


function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVDG})
    pvdg = parse_dead_pvt_table(f, outer_data)
    for tab in pvdg
        swap_unit_system_axes!(tab, units, (:pressure, :gas_formation_volume_factor, :viscosity))
    end
    data["PVDG"] = pvdg
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVTO})
    pvto = parse_live_pvt_table(f, outer_data)
    for tab in pvto
        swap_unit_system_axes!(tab["data"], units, (:pressure, :liquid_formation_volume_factor, :viscosity))
        swap_unit_system!(tab["key"], units, :u_rs)
    end
    data["PVTO"] = pvto
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVTG})
    pvtg = parse_live_pvt_table(f, outer_data)
    for tab in pvtg
        swap_unit_system_axes!(tab["data"], units, (:u_rv, :gas_formation_volume_factor, :viscosity))
        swap_unit_system!(tab["key"], units, :pressure)
    end
    data["PVTG"] = pvtg
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVTW})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN, NaN, NaN]
    utypes = (:pressure, :liquid_formation_volume_factor, :compressibility, :viscosity, :compressibility)
    nreg = number_of_tables(outer_data, :pvt)
    out = []
    for i = 1:nreg
        t = parse_defaulted_line(rec, tdims)
        swap_unit_system_axes!(t, units, utypes)
        @assert all(isfinite, t) "PVTW cannot be defaulted, found defaulted record in region $i"
        push!(out, t)
    end
    data["PVTW"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVCDO})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN, NaN, NaN]
    utypes = (:pressure, :liquid_formation_volume_factor, :compressibility, :viscosity, :compressibility)
    nreg = number_of_tables(outer_data, :pvt)
    out = []
    for i = 1:nreg
        t = parse_defaulted_line(rec, tdims)
        swap_unit_system_axes!(t, units, utypes)
        @assert all(isfinite, t) "PVCDO cannot be defaulted, found defaulted record in region $i"
        push!(out, t)
    end
    data["PVCDO"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVDO})
    pvdo = parse_dead_pvt_table(f, outer_data)
    for tab in pvdo
        swap_unit_system_axes!(tab, units, (:pressure, :liquid_formation_volume_factor, :viscosity))
    end
    data["PVDO"] = pvdo
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ROCK})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN, NaN, NaN, NaN]
    utypes = [:pressure, :compressibility, :compressibility, :compressibility, :id, :id]
    out = []
    nreg = number_of_tables(outer_data, :pvt)
    for i = 1:nreg
        l = parse_defaulted_line(rec, tdims)
        swap_unit_system_axes!(l, units, utypes)
        push!(out, l)
    end
    data["ROCK"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:COMPS})
    rec = read_record(f)
    ncomp = only(parse_defaulted_line(rec, [0]))
    data["COMPS"] = ncomp
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:DENSITY})
    tdims = [NaN, NaN, NaN]
    nreg = number_of_tables(outer_data, :pvt)
    out = []
    for i = 1:nreg
        rec = read_record(f)
        t = parse_defaulted_line(rec, tdims)
        swap_unit_system!(t, units, :density)
        push!(out, t)
    end
    data["DENSITY"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Union{Val{:RSCONST}, Val{:RSCONSTT}})
    rec = read_record(f)
    # TODO: This is missing units.
    tdims = [NaN, NaN]
    parsed = parse_defaulted_line(rec, tdims, required_num = length(tdims), keyword = "RSCONST")
    parser_message(cfg, outer_data, "RSCONST", PARSER_PARTIAL_SUPPORT)
    data["RSCONST"] = parsed
end

