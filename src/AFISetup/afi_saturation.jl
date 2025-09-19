function setup_saturation_variables(d::AFIInputFile, reservoir, sys; regions = setup_region_map(d))
    phase_symb = phases_symbol(sys)
    present = phases_present(sys)
    satfuns = find_records(d, "SaturationFunction", "IX", steps = false, model = true)
    krw = Dict()
    krg = Dict()
    krog = Dict()
    krow = Dict()
    krg = Dict()
    for satfun in map(x -> x, satfuns)
        @info "???" satfun.value satfun.keyword phase_symb
        vals = satfun.value
        reg = satfun.value["region"]
        if present.water
            krow[reg] = get_relperm(vals, "WaterRelPermFunction")
            if present.oil
                pc = get_saturation_function(vals, "CapPressure", "OilWaterCapPressureFunction")
                krow[reg] = get_relperm(vals, "OilInWaterRelPermFunction", reverse = false)
            end
        end

    end
end

function get_relperm(satfun, label; throw = true, reverse = false)
    krtab = get_saturation_function(satfun, "RelPerm", label, throw = throw)
    kr = krtab["RelPerm"]
    s = krtab["Saturation"]
    if reverse
        s = 1.0 .- s
    end
    return PhaseRelativePermeability(s, kr)
end

function get_saturation_function(satfun, type, label; throw = true)
    out = missing
    if haskey(satfun, label)
        table_name = satfun[label]
        return satfun[type][table_name]
    end
    if ismissing(out) && throw
        error("No $label for $type found in saturation function table.")
    end
    return out
end

function phases_present(sys)
    phases = JutulDarcy.get_phases(sys)
    has_water = AqueousPhase() in phases
    has_oil = LiquidPhase() in phases
    has_gas = VaporPhase() in phases
    return (water = has_water, oil = has_oil, gas = has_gas)
end

function phases_symbol(sys)
    present = phases_present(sys)
    has_water = present.water
    has_oil = present.oil
    has_gas = present.gas

    if has_water && has_oil && has_gas
        phases = :wog
    elseif has_water && has_oil
        phases = :wo
    elseif has_oil && has_gas
        phases = :og
    elseif has_water && has_gas
        phases = :wg
    elseif has_water
        phases = :w
    elseif has_oil
        phases = :o
    elseif has_gas
        phases = :g
    else
        error("No phases found in system")
    end
end