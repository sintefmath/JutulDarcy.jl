function setup_afi_schedule(afi::AFIInputFile, model::MultiModel)
    # well_constraints = Dict{Symbol, Any}()
    # well_status = Dict{Symbol, Any}()
    # well_type = Dict{Symbol, Any}()

    well_setup = Dict{Symbol, Any}()
    wells = get_model_wells(model)
    for (wname, w) in pairs(wells)
        well_setup[wname] = Dict{String, Any}(
            "Constraints" => Dict{String, Float64}(),
            "Status" => "OPEN",
            "Type" => "PRODUCER"
        )
    end

    reservoir = reservoir_domain(model)
    rmesh = physical_representation(reservoir)
    wells = get_model_wells(model)
    sys = reservoir_model(model).system
    # Constraints are in FM + completion status?
    # IX Steps:
    # WellDef:
    #  - connection_factor_reset
    #  - pi_multiplier_reset
    # WellToCellConnections status change

    ix_steps = afi["IX"]["STEPS"]
    fm_steps = afi["FM"]["STEPS"]
    dates = keys(ix_steps)
    @assert dates == keys(fm_steps)
    @assert issorted(dates)
    forces = []
    timesteps = Float64[]
    prev_date = missing
    for (dateno, date) in enumerate(dates)
        for (i, recv) in enumerate(ix_steps[date])
            kw = recv.keyword
            rec = recv.value
            # @info "IX $date: $i $kw" rec
            if kw == "WellDef"

            else
                println("IX: $kw: Skipping record $i for $date, record not used by setup code.")
            end
        end
        for (i, recv) in enumerate(fm_steps[date])
            kw = recv.keyword
            rec = recv.value
            # @info "FM $date: $i $kw" rec
            if kw == "Well"
                if haskey(rec, "name")
                    names = rec["name"]
                    if names isa AbstractString
                        names = [names]
                    end
                    for name in names
                        if contains(name, '*')
                            matcher = Regex(replace(name, '*' => ".*"))
                            for othername in keys(well_setup)
                                othername = String(othername)
                                if othername == name
                                    continue
                                end
                                is_match = !isnothing(findfirst(matcher, othername))
                                if is_match
                                    push!(names, othername)
                                end
                            end
                        end
                    end
                    new_names = String[]
                    for name in names
                        wname = Symbol(name)
                        wsetup = get(well_setup, wname, missing)
                        if ismissing(wsetup)
                            @warn "Well '$wname' in FM step $i for $date not found in model wells: $(keys(wells))"
                            continue
                        end
                        if haskey(rec, "Constraints")
                            verb = lowercase(rec["Constraints"].verb)
                            if verb != "add"
                                error("Only 'add' verb is implemented so far, got '$verb' in position 1 for record $i for $date: $(rec["Constraints"])")
                            end
                            for (cname, cval) in pairs(rec["Constraints"].constraints)
                                wsetup["Constraints"][cname] = Float64(cval)
                            end
                        end
                        if haskey(rec, "Status")
                            wsetup["Status"] = rec["Status"]
                        end
                        if haskey(rec, "Type")
                            wsetup["Type"] = rec["Type"]
                        end
                    end
                end
            else
                println("FM: $kw: Skipping record $i for $date, record not used by setup code.")
            end
        end
        if dateno > 1
            # Well constraints, etc
            f = forces_from_constraints(well_setup, date, sys, model, wells)
            push!(forces, f)
            # Timestep in seconds
            dt_i = date - prev_date
            dt_i::Millisecond
            push!(timesteps, Dates.value(dt_i)/1000.0)
        end
        prev_date = date
    end
    return (timesteps, forces)
end

function forces_from_constraints(well_setup, date, sys, model, wells)
    control = Dict{Symbol, Any}()
    limits = Dict{Symbol, Any}()
    phases = JutulDarcy.get_phases(sys)
    rhos = JutulDarcy.reference_densities(sys)

    function defaulted_constraint(c, default)
        x = get(c, default, missing)
        if ismissing(x)
            ks = keys(c)
            k = first(ks)
            x = c[k]
        else
            k = default
        end
        return (k, x)
    end
    for wname in keys(wells)
        wsetup = well_setup[wname]
        # @info "??" wname wsetup
        status = wsetup["Status"]
        wtype = lowercase(wsetup["Type"])
        c = wsetup["Constraints"]
        if status == GeoEnergyIO.IXParser.IX_OPEN
            if length(keys(c)) == 0
                println("No constraints found for well $wname at $date. Well will be disabled for step.")
                ctrl = JutulDarcy.setup_disabled_control()
            else
                if wtype == "producer"
                    ctrl_type_str, val = defaulted_constraint(c, "BOTTOM_HOLE_PRESSURE")
                    ctrl_type = control_type_to_symbol(ctrl_type_str)
                    ctrl = JutulDarcy.setup_producer_control(val, ctrl_type)
                    wsgn = -1.0
                else
                    if startswith(wtype, "water")
                        default = "WATER_INJECTION_RATE"
                        ph = AqueousPhase()
                    elseif startswith(wtype, "gas")
                        default = "GAS_INJECTION_RATE"
                        ph = VaporPhase()
                    elseif startswith(wtype, "oil")
                        default = "OIL_INJECTION_RATE"
                        ph = LiquidPhase()
                    else
                        error("Unknown injector type '$wtype' for well '$wname'")
                    end
                    ctrl_type_str, val = defaulted_constraint(c, default)
                    ctrl_type = control_type_to_symbol(ctrl_type_str)
                    ix = findfirst(isequal(ph), phases)
                    if isnothing(ix)
                        error("Phase '$ph' for injector well '$wname' not present in system phases: $phases")
                    end
                    mix = zeros(length(phases))
                    mix[ix] = 1.0
                    ctrl = JutulDarcy.setup_injector_control(val, ctrl_type, mix, density = rhos[ix])
                    wsgn = 1.0
                end
                lims = Dict()
                for (k, v) in pairs(c)
                    ck = control_type_to_symbol(k)
                    if endswith(uppercase(k), "_RATE")
                        v *= wsgn
                    end
                    lims[ck] = v
                end
                limits[wname] = lims
            end
        else
            ctrl = JutulDarcy.setup_disabled_control()
        end
        control[wname] = ctrl
    end
    return setup_reservoir_forces(model, control = control, limits = limits)
end

function control_type_to_symbol(s::String)
    s = lowercase(s)
    if s in ["bottom_hole_pressure", "bhp"]
        return :bhp
    elseif s in ["oil_production_rate", "oil_injection_rate", "orat"]
        return :orat
    elseif s in ["water_production_rate", "water_injection_rate", "wrat"]
        return :wrat
    elseif s in ["gas_production_rate", "gas_injection_rate", "grat"]
        return :grat
    elseif s in ["liquid_production_rate", "liquid_injection_rate", "lrat"]
        return :lrat
    else
        error("Unknown control type string '$s'")
    end
end
