function setup_afi_schedule(afi::AFIInputFile, model::MultiModel; step_limit = missing)
    OPEN = GeoEnergyIO.IXParser.IX_OPEN
    CLOSED = GeoEnergyIO.IXParser.IX_CLOSED
    well_setup = Dict{Symbol, Any}()
    group_lists = Dict{String, Any}(
        "StaticList" => Dict{String, Any}(),
        "Group" => Dict{String, Any}()
    )
    wells = get_model_wells(model)
    for (wname, w) in pairs(wells)
        orig_perf = model[wname].data_domain[:original_perforation_indices, JutulDarcy.Perforations()]
        remap = Dict{Int, Int}()
        for (i, p) in enumerate(orig_perf)
            remap[p] = i
        end
        wdomain = physical_representation(w)
        nperf = length(wdomain.perforations.reservoir)
        well_setup[wname] = Dict{String, Any}(
            "Constraints" => OrderedDict{String, Float64}(),
            "HistoryDataControl" => missing,
            "HistoricalControlModes" => missing,
            "Status" => "OPEN",
            "PiMultiplier" => ones(nperf),
            "ConnectionStatus" => fill(OPEN, nperf),
            "EffectivePiMultiplier" => ones(nperf),
            "Type" => "PRODUCER",
            "PerforationMap" => remap,
            "Enthalpy" => missing,
        )
    end

    stream_keys = ("FluidEnthalpy", "FluidStream", "FluidSourceExternal", "FluidSourceInternal")
    streams = Dict(k => Dict{String, Any}() for k in stream_keys)

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
    observation_data = obsh = get(afi["FM"], "OBSH", missing)
    dates = keys(ix_steps)
    @assert dates == keys(fm_steps)
    @assert issorted(dates)
    dates = collect(dates)
    forces = []
    timesteps = Float64[]
    t_elapsed = 0.0
    for (dateno, date) in enumerate(dates)
        for (i, recv) in enumerate(ix_steps[date])
            kw = recv.keyword
            rec = recv.value
            if kw == "WellDef"
                w2c = get(rec, "WellToCellConnections", missing)
                if !ismissing(w2c)
                    wname = Symbol(rec["WellName"])
                    pmap = well_setup[wname]["PerforationMap"]
                    eff_mult = well_setup[wname]["EffectivePiMultiplier"]
                    mult = well_setup[wname]["PiMultiplier"]
                    status = well_setup[wname]["ConnectionStatus"]

                    for (i, v) in enumerate(get(w2c, "PiMultiplier", []))
                        i_mapped = get(pmap, i, missing)
                        if ismissing(i_mapped)
                            continue
                        end
                        mult[i_mapped] = v
                    end
                    for (i, v) in enumerate(get(w2c, "Status", []))
                        i_mapped = get(pmap, i, missing)
                        if ismissing(i_mapped)
                            continue
                        end
                        status[pmap[i]] = v
                    end
                    for i in eachindex(eff_mult, mult, status)
                        next_mult = mult[i]
                        if status[i] != OPEN
                            next_mult = 0.0
                        end
                        eff_mult[i] = next_mult
                    end
                end
            else
                println("IX: $kw: Skipping record $i for $date, record not used by setup code.")
            end
        end
        for (i, recv) in enumerate(fm_steps[date])
            kw = recv.keyword
            rec = recv.value
            if kw == "StaticList" || kw == "Group"
                group_lists[kw][rec.group] = rec.members
            elseif kw == "HistoricalDataControl"
                if get(rec, "value", true)
                    wtype, wellnames = rec["Wells"]
                    @assert String(wtype) == "Well"
                    for wname in wellnames
                        wname = Symbol(wname)
                        well_setup[wname]["HistoryDataControl"] = rec["DataFileName"]
                    end
                end
            elseif kw == "Well"
                if haskey(rec, "name")
                    names = rec["name"]
                    if names isa AbstractString
                        names = [names]
                    end
                    new_names = String[]
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
                                    push!(new_names, othername)
                                end
                            end
                        else
                            push!(new_names, name)
                        end
                    end
                    names = new_names
                    for name in names
                        wname = Symbol(name)
                        wsetup = get(well_setup, wname, missing)
                        if ismissing(wsetup)
                            @warn "Well '$wname' in FM step $i for $date not found in model wells: $(keys(wells))"
                            continue
                        end
                        if haskey(rec, "Constraints")
                            rc = rec["Constraints"]
                            if rc isa GeoEnergyIO.IXParser.IXAssignmentRecord
                                cname = rc.index
                                cval = rc.value
                                wsetup["Constraints"][cname] = Float64(cval)
                            else
                                verb = lowercase(rc.verb)
                                if verb != "add"
                                    error("Only 'add' verb is implemented so far, got '$verb' in position 1 for record $i for $date: $(rec["Constraints"])")
                                end
                                for (cname, cval) in pairs(rc.constraints)
                                    wsetup["Constraints"][cname] = Float64(cval)
                                end
                            end
                        end
                        s = get(rec, "Status", "OPEN")
                        wsetup["Status"] = s
                        if s == "ALL_COMPLETIONS_SHUTIN" && false
                            well_setup[wname]["ConnectionStatus"] .= CLOSED
                            eff_mult = well_setup[wname]["EffectivePiMultiplier"]
                            eff_mult .= 0.0
                        end
                        if haskey(rec, "Type")
                            wsetup["Type"] = rec["Type"]
                        end
                        if haskey(rec, "Enthalpy")
                            wsetup["Enthalpy"] = rec["Enthalpy"].key
                        end
                        if haskey(rec, "HistoricalControlModes")
                            hcm = rec["HistoricalControlModes"]
                            if !(hcm isa AbstractVector)
                                hcm = [hcm]
                            end
                            wsetup["HistoricalControlModes"] = map(String, hcm)
                        end
                    end
                end
            elseif kw in stream_keys
                streams[kw][rec["group"]] = rec
            else
                println("FM: $kw: Skipping record $i for $date, record not used by setup code.")
            end
        end
        if dateno != length(dates)
            # Timestep in seconds
            dt_i = dates[dateno+1] - date
            dt_i::Millisecond
            dt_in_second = date_delta_to_seconds(dt_i)
            # Well constraints, etc
            f = forces_from_constraints(well_setup, observation_data, streams, date, sys, model, t_elapsed, dt_in_second, wells)
            push!(forces, f)
            push!(timesteps, dt_in_second)
            t_elapsed += dt_in_second
        end
        if !ismissing(step_limit) && dateno == step_limit
            break
        end
    end
    return (timesteps, forces)
end

function date_delta_to_seconds(dt_i::Millisecond)
    return Dates.value(dt_i)/1000.0
end

function forces_from_constraints(well_setup, observation_data, streams, date, sys, model, t_since_start_before_step, dt, wells)
    control = Dict{Symbol, Any}()
    limits = Dict{Symbol, Any}()
    phases = JutulDarcy.get_phases(sys)
    rhos = JutulDarcy.reference_densities(sys)
    rdomain = reservoir_domain(model)
    wforces = Dict{Symbol, Any}()

    for wname in keys(wells)
        wmodel = model[wname]
        wdomain = wmodel.data_domain
        wsetup = well_setup[wname]
        status = wsetup["Status"]
        wtype = lowercase(wsetup["Type"])
        c = wsetup["Constraints"]
        hctrl = wsetup["HistoryDataControl"]
        has_history = !ismissing(hctrl)
        if status == GeoEnergyIO.IXParser.IX_OPEN || status == "OPEN" || has_history
            if has_history
                hist_ctrl = wsetup["HistoricalControlModes"]
                ctrl_type, val, wtype = setup_history_control(hist_ctrl, wname, wtype, wsetup, observation_data, t_since_start_before_step, dt)
            else
                ctrl_type, val, wtype = setup_constraint_control(c, wtype)
            end
            if wtype == "disabled"
                ctrl = JutulDarcy.setup_disabled_control()
            elseif wtype == "producer"
                ctrl = JutulDarcy.setup_producer_control(val, ctrl_type)
                wsgn = -1.0
            else
                if startswith(wtype, "water")
                    ph = AqueousPhase()
                elseif startswith(wtype, "gas")
                    ph = VaporPhase()
                elseif startswith(wtype, "oil")
                    ph = LiquidPhase()
                else
                    error("Unknown injector type '$wtype' for well '$wname'")
                end
                ix = findfirst(isequal(ph), phases)
                if isnothing(ix)
                    error("Phase '$ph' for injector well '$wname' not present in system phases: $phases")
                end
                mix = zeros(length(phases))
                mix[ix] = 1.0
                enthalpy_stream = wsetup["Enthalpy"]
                if ismissing(enthalpy_stream)
                    T = convert_to_si(20.0, :Celsius)
                    enthalpy = missing
                else
                    enthalpy_info = streams["FluidEnthalpy"][enthalpy_stream]
                    T = convert_to_si(enthalpy_info["Temperature"], :Celsius)
                    # TODO: Handle regions here and treat enthalpy properly...
                    # P = enthalpy_info["Pressure"]
                    # rho_def = reservoir_model(model)[:PhaseMassDensities].pvt[ix]
                    # rho_c = rhos[ix]*JutulDarcy.shrinkage(rho_def, nothing, P, 1)
                    # wc = wmodel.domain.representation.perforations.reservoir
                    # hc = rdomain[:component_heat_capacity]
                    # if hc isa AbstractVector

                    # end
                    # heat_capacity = rdomain[:component_heat_capacity][ix, wc[1]]
                    # enthalpy = internal energy + p / rho
                    # internal energy = heat capacity * T
                    # enthalpy = (p, T) -> 
                    enthalpy = missing
                end
                ctrl = JutulDarcy.setup_injector_control(val, ctrl_type, mix,
                    density = rhos[ix],
                    temperature = T
                )
                wsgn = 1.0
            end
            if !(ctrl isa DisabledControl)
                lims = Dict()
                for (k, v) in pairs(c)
                    ck = control_type_to_symbol(k)
                    if ismissing(ck)
                        @warn "Unsupported constraint type '$k' for well '$wname', constraint will not be used."
                        continue
                    end
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
        mult = well_setup[wname]["EffectivePiMultiplier"]
        if any(x -> !(x ≈ 1.0), mult)
            mask = PerforationMask(mult)
            wforces[wname] = setup_forces(wmodel, mask = mask)
        end
        control[wname] = ctrl
    end
    f = setup_reservoir_forces(model; control = control, limits = limits, pairs(wforces)...)
    return f
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
    elseif s in ["injection_tubing_head_pressure"]
        return missing
    else
        error("Unknown control type string '$s'")
    end
end

function setup_constraint_control(c, wtype)
    function first_constraint(c)
        ks = keys(c)
        k = first(ks)
        x = c[k]
        return (k, x)
    end
    if length(c) == 0
        # println("No constraints provided for well, defaulting to disabled control.")
        wtype = "disabled"
        ctrl_type = :disabled
        val = missing
    else
        ctrl_type_str, val = first_constraint(c)
        ctrl_type = control_type_to_symbol(ctrl_type_str)
    end
    return (ctrl_type, val, wtype)
end

function setup_history_control(hist_ctrl, wname, wtype, wsetup, observation_data, t_since_start, dt)
    if ismissing(hist_ctrl)
        # println("HistoryDataControl is defaulted for well $wname but no constraint was provided as Well.HistoricalControlModes. Shutting well.")
        wtype = "disabled"
        return (:disabled, missing, wtype)
    elseif hist_ctrl isa AbstractVector
        if length(hist_ctrl) > 1
            println("Multiple HistoricalControlModes for well $wname at $date: $hist_ctrl, taking the first")
        end
        hmode = hist_ctrl[1]
    else
        hmode = hist_ctrl
    end
    function integrate_rate(k)
        feval = obs[k]
        nsub = 100
        dt_i = dt / nsub
        integrated_rate = 0.0
        for i in 1:nsub
            t_i = t_since_start + (i-0.5)*dt_i
            integrated_rate += feval(t_i)*dt_i
        end
        return integrated_rate/dt
    end
    ismissing(observation_data) && error("Observation data is required for HistoricalControlModes")
    maybe_ctrl(x, T) = ifelse(x ≈ 0.0, DisabledControl(), T(x))
    obs = observation_data[wsetup["HistoryDataControl"]]["wells_interp"][wname]
    if hmode == "BOTTOM_HOLE_PRESSURE"
        val = obs["BOTTOM_HOLE_PRESSURE"](t_since_start)
        ctrl_type = :bhp
    else
        if hmode == "RES_VOLUME_INJECTION_RATE"
            # println("RES_VOLUME_INJECTION_RATE handling not implemented yet, switching to rate constraint.")
            if wtype == "water_injector"
                ctrl_type, val, wtype = setup_history_control("WATER_INJECTION_RATE", wname, wtype, wsetup, observation_data, t_since_start, dt)
            elseif wtype == "gas_injector"
                ctrl_type, val, wtype = setup_history_control("GAS_INJECTION_RATE", wname, wtype, wsetup, observation_data, t_since_start, dt)
            else
                error("Unsupported well type '$wtype' for RES_VOLUME_INJECTION_RATE historical control for well '$wname'")
            end
        elseif hmode == "GAS_INJECTION_RATE"
            ctrl_type = :grat
            val = integrate_rate("GAS_INJECTION_RATE")
            wtype = "gas_injector"
        elseif hmode == "WATER_INJECTION_RATE"
            ctrl_type = :wrat
            val = integrate_rate("WATER_INJECTION_RATE")
            wtype = "water_injector"
        elseif hmode == "LIQUID_PRODUCTION_RATE"
            ctrl_type = :lrat
            val = -integrate_rate("LIQUID_PRODUCTION_RATE")
        elseif hmode == "GAS_PRODUCTION_RATE"
            val = -integrate_rate("GAS_PRODUCTION_RATE")
            ctrl_type = :grat
        elseif hmode == "RES_VOLUME_PRODUCTION_RATE"
            # println("RES_VOLUME_PRODUCTION_RATE handling not implemented yet, switching to LIQUID_PRODUCTION_RATE.")
            ctrl_type, val, wtype = setup_history_control("LIQUID_PRODUCTION_RATE", wname, wtype, wsetup, observation_data, t_since_start, dt)
        elseif hmode == "LIQUID_PRODUCTION_RATE"
            ctrl_type = :lrat
            val = -integrate_rate("LIQUID_PRODUCTION_RATE")
        elseif hmode == "OIL_PRODUCTION_RATE"
            ctrl_type = :orat
            val = -integrate_rate("OIL_PRODUCTION_RATE")
        else
            error("Unsupported HistoricalControlModes '$hmode' for well '$wname'")
        end
        if val ≈ 0.0
            wtype = "disabled"
        end
    end
    return (ctrl_type, val, wtype)
end
