function setup_afi_schedule(afi::AFIInputFile, model::MultiModel;
        step_limit = missing,
        step_override = missing,
        step_override_is_normalized = false,
        evaluate_enthalpy = true,
        injector_max_bhp_limit = Inf,
        producer_min_bhp_limit = si_unit(:atm)
    )
    CLOSED = GeoEnergyIO.IXParser.IX_CLOSED
    well_setup = Dict{Symbol, Any}()
    group_lists = Dict{String, Any}(
        "StaticList" => Dict{String, Any}(),
        "Group" => Dict{String, Any}()
    )
    wells = get_model_wells(model)
    for (wname, w) in pairs(wells)
        w_dd = model[wname].data_domain
        ids = w_dd[:perforation_connection_ids, JutulDarcy.Perforations()]
        nperf = length(physical_representation(w).perforations.reservoir)
        length(ids) == nperf && allunique(ids) ||
            throw(ArgumentError("Well $wname must have one unique AFI connection ID per perforation"))
        well_setup[wname] = Dict{String, Any}(
            "Constraints" => OrderedDict{String, Float64}(),
            "HistoryDataControl" => missing,
            "HistoricalControlModes" => missing,
            "Status" => "OPEN",
            "PiMultiplier" => ones(nperf),
            "ConnectionStatus" => fill(CLOSED, nperf),
            "Type" => "PRODUCER",
            "PerforationMap" => Dict(id => i for (i, id) in enumerate(ids)),
            "InputConnectionOrder" => OrderedDict{Tuple, Int}(),
            "Enthalpy" => missing,
            "injector_max_bhp_limit" => injector_max_bhp_limit,
            "producer_min_bhp_limit" => producer_min_bhp_limit
        )
    end

    completion_names = afi_completion_names(find_records(afi, "WellDef", "IX"))
    for record in find_records(afi, "WellDef", "IX", steps = false)
        update_afi_well_connections!(well_setup, record.value, completion_names)
    end

    stream_keys = ("FluidEnthalpy", "FluidStream", "FluidSourceExternal", "FluidSourceInternal")
    streams = Dict(k => Dict{String, Any}() for k in stream_keys)

    sys = reservoir_model(model).system
    is_geothermal = JutulDarcy.get_phases(sys) == (AqueousPhase(),) && JutulDarcy.model_is_thermal(model)
    if is_geothermal
        tab = JutulDarcy.Geothermal.geothermal_setup_tables()
    else
        tab = missing
    end
    ix_steps = afi["IX"]["STEPS"]
    fm_steps = afi["FM"]["STEPS"]
    observation_data = get(afi["FM"], "OBSH", missing)
    dates = keys(ix_steps)
    @assert dates == keys(fm_steps)
    @assert issorted(dates)
    dates = collect(dates)
    has_step_override = !ismissing(step_override)
    if has_step_override
        t_tot = date_delta_to_seconds(dates[end] - dates[1])
        if step_override_is_normalized
            step_override = step_override .* t_tot
        end
        if maximum(step_override) > t_tot + 1e-6
            @warn "Step override maximum time $(maximum(step_override)) exceeds total simulation time $t_tot"
        end
        if step_override[1] ≈ 0.0
            step_override = step_override[2:end]
        end
        all(diff(step_override) .> 0) || error("Step override times must be strictly increasing.")
    end
    forces = []
    timesteps = Float64[]
    t_elapsed = 0.0
    for (dateno, date) in enumerate(dates)
        for (i, recv) in enumerate(ix_steps[date])
            kw = recv.keyword
            rec = recv.value
            if kw == "WellDef"
                update_afi_well_connections!(well_setup, rec, completion_names)
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
                        ws = get(well_setup, wname, missing)
                        if ismissing(ws)
                            continue
                        end
                        ws["HistoryDataControl"] = rec["DataFileName"]
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
                        wsetup["Status"] = get(rec, "Status", wsetup["Status"])
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
            f = forces_from_constraints(well_setup, observation_data, streams, date, sys, model, t_elapsed, dt_in_second, wells, tab; evaluate_enthalpy = evaluate_enthalpy)
            if has_step_override
                for (idx, t_override) in enumerate(step_override)
                    inside_upper_bnd = t_override <= t_elapsed + dt_in_second
                    inside_lower_bnd = t_override > t_elapsed
                    if inside_lower_bnd && inside_upper_bnd
                        if idx == 1
                            dt_override = t_override
                        else
                            dt_override = t_override - step_override[idx-1]
                        end
                        push!(forces, f)
                        push!(timesteps, dt_override)
                    elseif !inside_upper_bnd
                        break
                    end
                end
            else
                push!(forces, f)
                push!(timesteps, dt_in_second)
            end
            t_elapsed += dt_in_second
        end
        if !ismissing(step_limit) && dateno == step_limit
            break
        end
    end
    if haskey(model.models, :Facility)
        active_wells = Symbol[]
        for f in forces
            for (wname, ctrl) in pairs(f[:Facility].control)
                if !(ctrl isa DisabledControl)
                    push!(active_wells, wname)
                end
            end
        end
        active_wells = unique!(active_wells)
        inactive_wells = setdiff(keys(wells), active_wells)
        ninactive = length(inactive_wells)
        if ninactive > 0
            println("$ninactive inactive wells in schedule (no provided constraints):\n$(join(inactive_wells, ", "))")
        end
    else
        println("No Facility model was set up. Case may be missing wells.")
    end
    return (timesteps, forces)
end

function update_afi_well_connections!(well_setup, record, completion_names)
    name = record["WellName"]
    well = get(well_setup, Symbol(name), nothing)
    isnothing(well) && return
    connections = get(record, "WellToCellConnections", Dict())
    ids = afi_connection_ids(name, connections, completion_names)
    input_order = well["InputConnectionOrder"]
    if isnothing(ids)
        # Indexed updates use AFI input order, including rows for inactive grid cells.
        positions = collect(values(input_order))
    else
        model_positions = well["PerforationMap"]
        positions = map(ids) do id
            get!(input_order, id) do
                i = get(model_positions, id, 0)
                if i > 0
                    well["ConnectionStatus"][i] = GeoEnergyIO.IXParser.IX_OPEN
                end
                return i
            end
        end
    end
    # Time-dependent Transmissibility and reset functions are not supported here.
    for field in ("PiMultiplier", "Status")
        updates = get(connections, field, missing)
        ismissing(updates) && continue
        target = field == "Status" ? well["ConnectionStatus"] : well["PiMultiplier"]
        for (row, value) in enumerate(updates)
            ismissing(value) && continue
            row <= length(positions) ||
                throw(ArgumentError("Well $name has a $field update at row $row with no corresponding AFI connection"))
            i = positions[row]
            i == 0 && continue
            target[i] = field == "PiMultiplier" ? max(value, 0.0) : value
        end
    end
    return nothing
end

function date_delta_to_seconds(dt_i::Millisecond)
    return Dates.value(dt_i)/1000.0
end

function forces_from_constraints(well_setup, observation_data, streams, date, sys, model, t_since_start_before_step, dt, wells, tab; evaluate_enthalpy = false)
    control = Dict{Symbol, Any}()
    limits = Dict{Symbol, Any}()
    phases = JutulDarcy.get_phases(sys)
    rhos = JutulDarcy.reference_densities(sys)
    rdomain = reservoir_domain(model)
    wforces = Dict{Symbol, Any}()

    for wname in keys(wells)
        wmodel = model[wname]
        wsetup = well_setup[wname]
        status = wsetup["Status"]
        open = wsetup["ConnectionStatus"] .== GeoEnergyIO.IXParser.IX_OPEN
        mult = ifelse.(open, wsetup["PiMultiplier"], 0.0)
        if status == "ALL_COMPLETIONS_SHUTIN"
            fill!(mult, 0.0)
        end
        if any(x -> !(x ≈ 1.0), mult)
            wforces[wname] = setup_forces(wmodel, mask = PerforationMask(mult))
        end
        do_shut = all(iszero, mult)
        wtype = lowercase(wsetup["Type"])
        c = wsetup["Constraints"]
        hctrl = wsetup["HistoryDataControl"]
        has_history = !ismissing(hctrl)
        if status == GeoEnergyIO.IXParser.IX_OPEN || status == "OPEN" || has_history
            if has_history
                hist_ctrl = wsetup["HistoricalControlModes"]
                ctrl_type, val, wtype = setup_history_control(hist_ctrl, wname, wtype, wsetup, observation_data, t_since_start_before_step, c, dt)
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
                    if ismissing(tab) || !evaluate_enthalpy
                        enthalpy = missing
                    else
                        p = enthalpy_info["Pressure"]
                        # Consistent with the tabulated enthalpy used in
                        # the model
                        enthalpy = tab[:enthalpy](p, T)
                    end
                end
                ctrl = JutulDarcy.setup_injector_control(val, ctrl_type, mix,
                    density = rhos[ix],
                    temperature = T,
                    enthalpy = enthalpy
                )
                wsgn = 1.0
            end
            if do_shut && !(ctrl isa DisabledControl)
                println("Well $wname has all perforations closed at $date, disabling well for this step.")
                ctrl = JutulDarcy.setup_disabled_control()
            end
            control[wname] = ctrl
            ibhp_max = wsetup["injector_max_bhp_limit"]
            pbhp_min = wsetup["producer_min_bhp_limit"]
            if !(ctrl isa DisabledControl)
                is_injector = ctrl isa InjectorControl
                lims = Dict()
                for (k, v) in pairs(c)
                    ck = control_type_to_symbol(k)
                    if ismissing(ck)
                        @warn "Unsupported constraint type '$k' for well '$wname', constraint will not be used."
                        continue
                    end
                    if endswith(uppercase(k), "_RATE")
                        if abs(v) < 1e-16
                            # Do not set zero limits!
                            continue
                        end
                        v *= wsgn
                    end
                    if ck == :bhp
                        if is_injector
                            v = min(v, ibhp_max)
                        else
                            v = max(v, pbhp_min)
                        end
                    end
                    lims[ck] = v
                end
                limits[wname] = lims
            end
        else
            control[wname] = JutulDarcy.setup_disabled_control()
        end
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

function setup_history_control(hist_ctrl, wname, wtype, wsetup, observation_data, t_since_start, constraints, dt)
    if ismissing(hist_ctrl)
        if wtype == "producer"
            hmode = "LIQUID_PRODUCTION_RATE"
        elseif wtype == "water_injector"
            hmode = "WATER_INJECTION_RATE"
        elseif wtype == "gas_injector"
            hmode = "GAS_INJECTION_RATE"
        else
            error("Unsupported well type '$wtype' for historical control for well '$wname'")
        end
        println("$wname ($wtype): Under HistoryDataControl, no constraint in Well.HistoricalControlModes. Setting to $hmode.")
    elseif hist_ctrl isa AbstractVector
        if length(hist_ctrl) > 1
            println("Multiple HistoricalControlModes for well $wname at $t_since_start s: $hist_ctrl, taking the first")
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
    function get_constraint_value(k)
        # Constraints override the integrated rate from observation data
        cval = get(constraints, k, 0.0)
        if cval == 0.0
            cval = integrate_rate(k)
        end
        return cval
    end
    ismissing(observation_data) && error("Observation data is required for HistoricalControlModes")
    maybe_ctrl(x, T) = ifelse(x ≈ 0.0, DisabledControl(), T(x))
    obs = get(observation_data[wsetup["HistoryDataControl"]]["wells_interp"], wname, missing)
    if ismissing(obs)
        jutul_message("OBSH", "No observation data found for well '$wname' with HistoricalDataControl", color = :red)
        ctrl_type = :disabled
        val = 0.0
        wtype = "disabled"
    elseif hmode == "BOTTOM_HOLE_PRESSURE"
        val = obs["BOTTOM_HOLE_PRESSURE"](t_since_start)
        ctrl_type = :bhp
    else
        if hmode == "RES_VOLUME_INJECTION_RATE"
            # println("RES_VOLUME_INJECTION_RATE handling not implemented yet, switching to rate constraint.")
            if wtype == "water_injector"
                ctrl_type, val, wtype = setup_history_control("WATER_INJECTION_RATE", wname, wtype, wsetup, observation_data, t_since_start, constraints, dt)
            elseif wtype == "gas_injector"
                ctrl_type, val, wtype = setup_history_control("GAS_INJECTION_RATE", wname, wtype, wsetup, observation_data, t_since_start, constraints, dt)
            else
                error("Unsupported well type '$wtype' for RES_VOLUME_INJECTION_RATE historical control for well '$wname'")
            end
        elseif hmode == "GAS_INJECTION_RATE"
            ctrl_type = :grat
            val = get_constraint_value("GAS_INJECTION_RATE")
            wtype = "gas_injector"
        elseif hmode == "WATER_INJECTION_RATE"
            ctrl_type = :wrat
            val = get_constraint_value("WATER_INJECTION_RATE")
            wtype = "water_injector"
        elseif hmode == "LIQUID_PRODUCTION_RATE"
            ctrl_type = :lrat
            val = -get_constraint_value("LIQUID_PRODUCTION_RATE")
        elseif hmode == "GAS_PRODUCTION_RATE"
            val = -get_constraint_value("GAS_PRODUCTION_RATE")
            ctrl_type = :grat
        elseif hmode == "RES_VOLUME_PRODUCTION_RATE"
            ctrl_type, val, wtype = setup_history_control("LIQUID_PRODUCTION_RATE", wname, wtype, wsetup, observation_data, t_since_start, constraints, dt)
        elseif hmode == "LIQUID_PRODUCTION_RATE"
            ctrl_type = :lrat
            val = -get_constraint_value("LIQUID_PRODUCTION_RATE")
        elseif hmode == "OIL_PRODUCTION_RATE"
            ctrl_type = :orat
            val = -get_constraint_value("OIL_PRODUCTION_RATE")
        else
            error("Unsupported HistoricalControlModes '$hmode' for well '$wname'")
        end
        if val ≈ 0.0
            wtype = "disabled"
        end
    end
    return (ctrl_type, val, wtype)
end
