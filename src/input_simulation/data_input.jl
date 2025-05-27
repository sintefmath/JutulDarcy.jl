"""
    simulate_data_file(inp; parse_arg = NamedTuple(), kwarg...)

Simulate standard input file (with extension .DATA). `inp` can either be the
output from the GeoEnergyIO function `parse_data_file` or a `String` for the
path of an input file with the .DATA extension.

Additional arguments are passed onto [`simulate_reservoir`](@ref). Extra inputs
to the parser can be sent as a `setup_arg` `NamedTuple`.
"""
function simulate_data_file(data; setup_arg = missing, mode = :default, kwarg...)
    if ismissing(setup_arg)
        setup_arg = (split_wells = mode != :default, )
    end
    case = setup_case_from_data_file(data; setup_arg...)
    result = simulate_reservoir(case; mode = mode, kwarg...)
    result.extra[:case] = case
    return result
end

"""
    case = setup_case_from_data_file(
        filename;
        parse_arg = NamedTuple(),
        kwarg...
    )

    case, data = setup_case_from_data_file(filename; include_data = true)

Set up a [`JutulCase`](@ref) from a standard input file (with extension .DATA).
Additional arguments to the parser can be passed as key/values in a `NamedTuple`
given as `parse_arg`. The optional input `include_data=true` will make the
function return the parsed data in addition the case itself. Additional keyword
arguments are passed on to [`setup_case_from_parsed_data`](@ref).
"""
function setup_case_from_data_file(
        fn;
        parse_arg = NamedTuple(),
        include_data::Bool = false,
        kwarg...
    )
    if fn isa String
        data = parse_data_file(fn; parse_arg...)
    else
        @assert fn isa AbstractDict
        data = fn
    end
    case = setup_case_from_parsed_data(data; kwarg...)
    if include_data
        out = (case, data)
    else
        out = case
    end
    return out
end

"""
    setup_case_from_parsed_data(datafile; skip_wells = false, simple_well = true, use_ijk_trans = true, verbose = false, kwarg...)

Set up a case from a parsed input file (in `Dict` format). Internal function,
not exported. Use [`setup_case_from_data_file`](@ref).
"""
function setup_case_from_parsed_data(datafile;
        skip_wells = false,
        skip_forces = skip_wells,
        simple_well = true,
        use_ijk_trans = true,
        verbose = false,
        zcorn_depths = true,
        normalize = true,
        convert_co2store = true,
        repair_zcorn = true,
        process_pinch = true,
        kwarg...
    )
    function msg(s)
        if verbose
            jutul_message("Setup", s)
        end
    end
    msg("Parsing physics and system.")
    rs = datafile["RUNSPEC"]
    can_convert_co2 = haskey(rs, "CO2STORE") || haskey(rs, "CO2STOR")
    if convert_co2store && can_convert_co2
        msg("CO2STORE/CO2STOR found, calling converter...")
        datafile = convert_co2store_to_co2_brine(datafile)
        rs = datafile["RUNSPEC"]
    end
    sys, pvt = parse_physics_types(datafile, pvt_region = 1)
    is_blackoil = sys isa StandardBlackOilSystem
    is_compositional = sys isa CompositionalSystem || sys == :co2brine
    is_thermal = haskey(datafile["RUNSPEC"], "THERMAL")

    water = haskey(rs, "WATER")
    if is_compositional
        oil = gas = true
    else
        oil = haskey(rs, "OIL")
        gas = haskey(rs, "GAS")
    end
    if haskey(datafile, "PROPS")
        props = datafile["PROPS"]
    else
        props = Dict{String, Any}()
    end

    msg("Parsing reservoir domain.")
    domain = parse_reservoir(datafile, zcorn_depths = zcorn_depths, repair_zcorn = repair_zcorn, process_pinch = process_pinch)
    pvt_reg = reservoir_regions(domain, :pvtnum)
    has_pvt = isnothing(pvt_reg)
    # Parse wells
    msg("Parsing schedule.")
    if skip_forces
        wells = []
        controls = limits = cstep = well_forces = missing
        dt = [1.0]
    else
        wells, controls, limits, cstep, dt, well_forces = parse_schedule(domain, sys, datafile; simple_well = simple_well)
        if skip_wells
            empty!(wells)
        end
    end
    msg("Setting up model with $(length(wells)) wells.")
    wells_pvt = Dict()
    wells_systems = []
    for w in wells
        if has_pvt
            c = first(w.perforations.reservoir)
            reg_w = pvt_reg[c]
            sys, pvt = parse_physics_types(datafile, pvt_region = reg_w)
        else
            sys_w, pvt_w = sys, pvt
        end
        wells_pvt[w.name] = pvt_w
        push!(wells_systems, sys_w)
    end

    extra_arg = Dict{Symbol, Any}()
    if haskey(rs, "SALTS")
        extra_arg[:salt_mole_fractions] = rs["SALT_MOLE_FRACTIONS"]
        extra_arg[:salt_names] = rs["SALTS"]
    end
    tracers = data_file_active_tracers(datafile)
    if length(tracers) > 0
        tracer_types = []
        for t in tracers
            if t == "POLYMER"
                push!(tracer_types, Tracers.PolymerTracer(sys))
                # push!(tracer_types, MultiPhaseTracer(sys))
            else
                error("Tracer $t not supported.")
            end
        end
        extra_arg[:tracers] = [t for t in tracer_types]
    end

    model = setup_reservoir_model(domain, sys;
        thermal = is_thermal,
        wells = wells,
        extra_out = false,
        wells_systems = wells_systems,
        kwarg...,
        extra_arg...
    )
    for (k, submodel) in pairs(model.models)
        if model_or_domain_is_well(submodel) || k == :Reservoir
            # Modify secondary variables
            if !is_compositional
                svar = submodel.secondary_variables
                pvt_reg_i = reservoir_regions(submodel.data_domain, :pvtnum)
                if model_or_domain_is_well(submodel)
                    pvt_i = wells_pvt[submodel.domain.representation.name]
                else
                    pvt_i = pvt
                end
                pvt_i = tuple(pvt_i...)

                if is_thermal && haskey(props, "WATDENT")
                    watdent = WATDENT(props["WATDENT"])
                else
                    watdent = nothing
                end
                rho = DeckPhaseMassDensities(pvt_i, regions = pvt_reg_i, watdent = watdent)
                if sys isa StandardBlackOilSystem
                    b_i = DeckShrinkageFactors(pvt_i, regions = pvt_reg_i, watdent = watdent)
                    set_secondary_variables!(submodel,
                        ShrinkageFactors = wrap_reservoir_variable(sys, b_i, :flow)
                    )
                end
                if is_thermal
                    if haskey(props, "WATVISCT") || haskey(props, "OILVISCT") || haskey(props, "GASVISCT")
                        thermal_visc = DeckThermalViscosityTable(props, pvt_i, water, oil, gas)
                    else
                        thermal_visc = nothing
                    end
                else
                    thermal_visc = nothing
                end
                mu = DeckPhaseViscosities(pvt_i, regions = pvt_reg_i, thermal = thermal_visc)
                set_secondary_variables!(submodel,
                    PhaseViscosities = wrap_reservoir_variable(sys, mu, :flow),
                    PhaseMassDensities = wrap_reservoir_variable(sys, rho, :flow)
                )
            end
            if is_thermal
                set_thermal_deck_specialization!(submodel, props, domain[:pvtnum], oil, water, gas)
            end
            if k == :Reservoir
                set_deck_specialization!(submodel, rs, props, domain[:satnum], oil, water, gas)
            end
        end
    end
    if "POLYMER" in tracers
        Tracers.set_polymer_model!(model, datafile)
    end
    msg("Setting up forces.")
    if skip_forces
        forces = setup_reservoir_forces(model)
    else
        forces = parse_forces(model, datafile, sys, wells, controls, limits, cstep, dt, well_forces)
    end
    msg("Setting up initial state.")
    state0 = parse_state0(model, datafile, normalize = normalize)
    msg("Setting up parameters.")
    parameters = setup_parameters(model)
    if haskey(props, "SWL")
        G = physical_representation(domain)
        swl = vec(props["SWL"])
        # parameters[:Reservoir][:ConnateWater] .= swl[G.cell_map]
    end
    if use_ijk_trans
        parameters[:Reservoir][:Transmissibilities] = reservoir_transmissibility(domain, version = :ijk);
    end
    msg("Setup complete.")
    return JutulCase(model, dt, forces, state0 = state0, parameters = parameters, input_data = datafile)
end

function parse_well_from_compdat(domain, wname, cdat, wspecs, msdata, compord, step; simple_well = isnothing(msdata))
    wc, WI, open, = compdat_to_connection_factors(domain, wspecs, cdat, step, sort = true, order = compord)
    ref_depth = wspecs.ref_depth
    if isnan(ref_depth)
        ref_depth = nothing
    end
    W = missing
    accumulator_volume = missing
    if simple_well
        @assert isnothing(msdata)
        accumulator_volume = 0.05*mean(domain[:volumes][wc])*mean(domain[:porosity][wc])
    else
        if !isnothing(msdata)
            has_welsegs = haskey(msdata, "WELSEGS")
            has_compsegs = haskey(msdata, "COMPSEGS")
            if has_welsegs || has_compsegs
                @assert has_welsegs && has_compsegs "Both COMPSEGS and WELSEGS must be defined"
                @assert msdata["COMPSEGS"][1] == wname
                @assert msdata["WELSEGS"][1] == wname

                welsegs = msdata["WELSEGS"][2]
                header = welsegs.header
                segments = welsegs.segments

                compsegs = msdata["COMPSEGS"][2]

                top_depth = header[2]
                if isnothing(ref_depth)
                    ref_depth = top_depth
                elseif !(ref_depth ≈ top_depth)
                    @warn "Reference depths for ms well should coincide with top depth: ref depth $ref_depth != top depth $top_depth"
                end
                top_tubing_delta = header[3]
                accumulator_volume = header[4]
                segment_increment_type = header[5]
                top_x = header[8]
                top_y = header[9]
                @assert segment_increment_type in ("ABS", "INC")
                is_inc = segment_increment_type == "INC"
                max_seg = maximum(x -> x[1], segments, init = 0)
                num_edges = max_seg-1
                conn = Tuple{Int, Int}[]

                top_node = [top_x, top_y, top_depth]
                centers = zeros(3, num_edges)
                volumes = zeros(num_edges)
                diameter = fill(NaN, num_edges)
                branches = zeros(Int, num_edges)
                tubing_lengths = zeros(num_edges)
                tubing_depths = zeros(num_edges)
                # tubing_depths[1] = top_tubing_delta

                segment_models = Vector{SegmentWellBoreFrictionHB{Float64}}(undef, num_edges)
                for segment in segments
                    start, stop, branch, start_conn,
                    dist, depth_delta, D, rough, cross_sect, vol, = segment
                    parts_in_segment = stop - start + 1
                    if start_conn == 1
                        depth_at_start = top_depth
                        tubing_at_start = top_tubing_delta
                    else
                        depth_at_start = centers[3, start_conn-1]
                        tubing_at_start = tubing_depths[start_conn-1]
                    end
                    if is_inc
                        L = dist
                        Δz = depth_delta
                    else
                        L = (dist - tubing_at_start)/parts_in_segment
                        Δz = (depth_delta - depth_at_start)/parts_in_segment
                    end
                    if cross_sect <= 0.0
                        cross_sect = π*(D/2)^2
                    end
                    if vol <= 0.0
                        vol = cross_sect*L
                    end
                    Δp = SegmentWellBoreFrictionHB(L, rough, D)

                    push!(conn, (start, start_conn))
                    current_tubing = tubing_at_start
                    current_depth = depth_at_start
                    for (ix, seg_ix) in enumerate(start:stop)
                        current_tubing += L
                        current_depth += Δz

                        edge_ix = seg_ix - 1
                        @assert isnan(diameter[edge_ix]) "Values are being overwritten in ms well - programming error?"
                        segment_models[edge_ix] = Δp
                        diameter[edge_ix] = D
                        volumes[edge_ix] = vol
                        branches[edge_ix] = branch
                        tubing_lengths[edge_ix] = L
                        # TODO: Set better x, y
                        centers[1, edge_ix] = top_x
                        centers[2, edge_ix] = top_y
                        centers[3, edge_ix] = current_depth
                        tubing_depths[edge_ix] = current_tubing
                    end
                end
                N = zeros(Int, 2, length(conn))
                for (i, t) in enumerate(conn)
                    N[:, i] .= t
                end
                perforation_cells = map_compdat_to_multisegment_segments(compsegs, branches, tubing_lengths, tubing_depths, cdat)
                cell_centers = domain[:cell_centroids]
                dz = vec(cell_centers[3, wc]) - vec(centers[3, perforation_cells])
                W = MultiSegmentWell(wc, volumes, centers;
                    name = Symbol(wname),
                    WI = WI,
                    dz = dz, # TODO
                    N = N,
                    perforation_cells = perforation_cells,
                    reference_depth = ref_depth,
                    segment_models = segment_models,
                )
            end
        end
    end

    if ismissing(W)
        W = setup_well(domain, wc;
            name = Symbol(wname),
            accumulator_volume = accumulator_volume,
            WI = WI,
            reference_depth = ref_depth,
            simple_well = simple_well
        )
    end
    return (W, wc, WI, open)
end

function map_compdat_to_multisegment_segments(compsegs, branches, tubing_lengths, tubing_depths, cdat)
    completions = keys(cdat)
    perforation_cells = zeros(Int, length(completions))
    segment_ijk = map(x -> (x[1], x[2], x[3]), compsegs)
    for (completion_index, completion) in enumerate(completions)
        segment_candidates = findall(isequal(completion), segment_ijk)
        if length(segment_candidates) == 0
            @warn "No segments found for completion $completion. Expanding search to all segments."
            segment_candidates = eachindex(compsegs)
        end
        prev_dist = Inf
        closest = -1
        for segment_index in segment_candidates
            I, J, K, branch, tube_start, tube_end, dir, dir_ijk, cdepth, = compsegs[segment_index]
            for (i, b) in enumerate(branches)
                if b == branch
                    L = tubing_lengths[i]
                    seg_end = tubing_depths[i]
                    seg_mid = seg_end - L/2
                    tube_mid = (tube_end + tube_start)/2
                    dist = abs(seg_mid - tube_mid)
                    if dist < prev_dist
                        closest = i
                        prev_dist = dist
                    end
                end
            end
        end
        perforation_cells[completion_index] = closest
    end
    return perforation_cells
end

function compdat_to_connection_factors(domain, wspec, v, step; sort = true, order = "TRACK", ijk_lookup = missing)
    G = physical_representation(domain)
    if ismissing(ijk_lookup)
        ijk_lookup = ijk_lookup_dict(G)
    end
    K = domain[:permeability]
    T = eltype(K)
    if haskey(domain, :net_to_gross)
        T = promote_type(T, eltype(domain[:net_to_gross]))
    end

    ij_ix = collect(keys(v))
    wc = map(i -> ijk_lookup[i], ij_ix)

    getf(k) = map(x -> v[x][k], ij_ix)
    open = getf(:open)
    d = T.(getf(:diameter))
    Kh = T.(getf(:Kh))
    WI = T.(getf(:WI))
    skin = T.(getf(:skin))
    dir = getf(:dir)
    mul = getf(:mul)
    fresh = map(x -> v[x].ctrl == step, ij_ix)

    for i in eachindex(WI)
        W_i = WI[i]
        if isnan(W_i) || W_i < 0.0
            c = wc[i]
            if K isa AbstractVector
                k_i = K[c]
            else
                k_i = K[:, c]
            end
            drainage_radius = wspec.drainage_radius
            if drainage_radius <= 0.0
                drainage_radius = nothing
            end
            if haskey(domain, :net_to_gross)
                ntg = domain[:net_to_gross][c]
            else
                ntg = 1.0
            end
            WI[i] = compute_peaceman_index(G, k_i, d[i]/2, c,
                skin = skin[i],
                Kh = Kh[i],
                net_to_gross = ntg,
                dir = Symbol(dir[i]),
                drainage_radius = drainage_radius
            )
        end
    end
    if sort
        ix = well_completion_sortperm(domain, wspec, order, wc, dir)
        wc = wc[ix]
        WI = WI[ix]
        open = open[ix]
        fresh = fresh[ix]
    end
    return (wc, WI, open, mul, fresh)
end

function parse_schedule(domain, runspec, props, schedule, sys; simple_well = true)
    G = physical_representation(domain)
    ijk_lookup = ijk_lookup_dict(G)
    dt, cstep, controls, status, completions, msdata, limits = parse_control_steps(runspec, props, schedule, sys)
    completions, bad_wells = filter_inactive_completions!(completions, G, ijk_lookup)
    @assert length(controls) == length(completions)
    handle_wells_without_active_perforations!(bad_wells, completions, controls, limits)
    ncomp = length(completions)
    wells = []
    well_forces = Dict{Symbol, Any}[]
    for i in eachindex(completions)
        push!(well_forces, Dict{Symbol, Any}())
    end
    for (k, v) in pairs(completions[end])
        wspec = schedule["WELSPECS"][k]
        compord = schedule["COMPORD"][k]
        for i in eachindex(completions)
            well_forces[i][Symbol(k)] = (mask = nothing, )
        end
        if haskey(msdata, k)
            msdata_k = msdata[k]
            k_is_simple_well = false
        else
            msdata_k = nothing
            k_is_simple_well = simple_well
        end
        W, wc_base, WI_base, open = parse_well_from_compdat(domain, k, v, wspec, msdata_k, compord, length(completions); simple_well = k_is_simple_well)
        n_wi = length(WI_base)
        T = eltype(WI_base)
        wpi_mul = ones(T, n_wi)
        for (i, c) in enumerate(completions)
            compdat = c[k]
            well_control = controls[i][k]
            flag = status[i][k]
            well_is_shut = well_control isa DisabledControl
            wi_mul = zeros(T, n_wi)
            if !well_is_shut
                wc, WI, open, multipliers, fresh = compdat_to_connection_factors(domain, wspec, compdat, i, ijk_lookup = ijk_lookup, sort = false)
                for (c, wi, is_open, is_fresh, mul) in zip(wc, WI, open, fresh, multipliers)
                    compl_idx = findfirst(isequal(c), wc_base)
                    if isnothing(compl_idx)
                        # This perforation is missing, leave at zero.
                        continue
                    end
                    if !is_fresh
                        # Perforation is not fresh and could have previously
                        # gotten a WPIMULT applied. TODO: Some of the logic
                        # around WPIMULT is not 100% clear. The current
                        # implementation should be fine up to the stage where
                        # many different keywords start to interact for the same
                        # timestep.
                        mul *= wpi_mul[compl_idx]
                    end
                    wpi_mul[compl_idx] = mul
                    wi_mul[compl_idx] = wpi_mul[compl_idx]*wi*is_open/WI_base[compl_idx]
                end
            end
            if !well_is_shut
                # Multiply by well factor to account for downtime in
                # pressure drop when well is not always active
                @. wi_mul *= well_control.factor
            end
            if flag == "SHUT"
                @. wi_mul = 0.0
            end
            if all(isapprox(1.0), wi_mul)
                mask = nothing
            else
                mask = PerforationMask(wi_mul)
            end
            # TODO: Make this use setup_forces properly
            well_forces[i][Symbol(k)] = (mask = mask, )
        end
        push!(wells, W)
    end
    return (wells, controls, limits, cstep, dt, well_forces)
end

function ijk_lookup_dict(G)
    ijk_lookup = Dict{Tuple{Int, Int, Int}, Int}()
    for i in 1:number_of_cells(G)
        ijk = cell_ijk(G, i)
        ijk_lookup[ijk] = i
    end
    return ijk_lookup
end

function filter_inactive_completions!(completions_vector, g, ijk_lookup)
    tuple_active = Dict{Tuple{String, NTuple{3, Int}}, Bool}()
    for completions in completions_vector
        for (well, completion_set) in pairs(completions)
            for tupl in keys(completion_set)
                k = (well, tupl)
                if !haskey(tuple_active, k)
                    if haskey(ijk_lookup, tupl)
                        tuple_active[k] = true
                    else
                        jutul_message("$well completion",
                            "Removed COMPDAT as $tupl is not active in processed mesh.",
                            color = :yellow
                        )
                        tuple_active[k] = false
                    end
                end
                if !tuple_active[k]
                    delete!(completion_set, tupl)
                end
            end
        end
    end
    bad_wells = String[]
    for (well, completion_set) in pairs(completions_vector[end])
        if length(keys(completion_set)) == 0
            for c in completions_vector
                @assert length(keys(c[well])) == 0
            end
            push!(bad_wells, well)
        end
    end
    return (completions_vector, bad_wells)
end

function handle_wells_without_active_perforations!(bad_wells, completions, controls, limits)
    if length(bad_wells) > 0
        @assert length(completions) == length(controls) == length(limits)
        for well in bad_wells
            println("$well has no completions in active cells. Well will be not be present in model or simulation results.")
            for i in eachindex(completions)
                delete!(completions[i], well)
                delete!(controls[i], well)
                delete!(limits[i], well)
            end
        end
    end
end

function parse_forces(model, datafile, sys, wells, controls, limits, cstep, dt, well_forces)
    bc, sources = parse_aquifer_bc(model, datafile, sys)
    if length(wells) == 0
        return setup_reservoir_forces(model, bc = bc, sources = sources)
    end
    forces = []
    @assert length(controls) == length(limits) == length(well_forces)
    i = 0
    for (ctrl, lim, wforce) in zip(controls, limits, well_forces)
        i += 1
        ctrl_s = Dict{Symbol, Any}()
        for (k, v) in pairs(ctrl)
            wf_k = wforce[Symbol(k)]
            if !isa(v, DisabledControl) && !isnothing(wf_k.mask)
                if all(isequal(0), wf_k.mask.values)
                    jutul_message("Shutting $k", "Well has no open perforations at step $i, shutting.")
                    v = DisabledControl()
                end
            end
            ctrl_s[Symbol(k)] = v
        end
        lim_s = Dict{Symbol, Any}()
        for (k, v) in pairs(lim)
            lim_s[Symbol(k)] = v
        end
        f = setup_reservoir_forces(model;
            bc = bc,
            sources = sources,
            control = ctrl_s,
            limits = lim_s,
            pairs(wforce)...)
        push!(forces, f)
    end
    return forces[cstep]
end

function parse_state0(model, datafile; normalize = true)
    rmodel = reservoir_model(model)
    reservoir = reservoir_domain(rmodel)
    init = Dict{Symbol, Any}()
    sol = datafile["SOLUTION"]

    if haskey(sol, "EQUIL")
        init = parse_state0_equil(rmodel, datafile; normalize = normalize)
    else
        init = parse_state0_direct_assignment(rmodel, datafile)
    end
    if haskey(reservoir, :numerical_aquifers)
        # Aquifers are a special case
        initialize_numerical_aquifers!(init, rmodel, reservoir[:numerical_aquifers])
    end

    if haskey(init, :Temperature)
        # Temperature can be set during equil, write it to data domain for
        # parameter initialization.
        rmodel.data_domain[:temperature] = init[:Temperature]
        if !haskey(rmodel.primary_variables, :Temperature)
            delete!(init, :Temperature)
        end
    end
    state0 = setup_reservoir_state(model, init)
    return state0
end

function parse_state0_direct_assignment(model, datafile)
    sys = model.system
    init = Dict{Symbol, Any}()
    sol = datafile["SOLUTION"]
    G = physical_representation(model.data_domain)
    nc = number_of_cells(G)
    ix = G.cell_map
    is_blackoil = sys isa StandardBlackOilSystem
    is_compositional = sys isa MultiPhaseCompositionalSystemLV

    function get_active(k; to_zero = false)
        if haskey(sol, k)
            x = zeros(nc)
            for (i, c) in enumerate(G.cell_map)
                x[i] = sol[k][c]
            end
        elseif to_zero
            x = zeros(nc)
        else
            x = missing
        end
        return x
    end

    function get_active_fraction(k)
        ϵ = MultiComponentFlash.MINIMUM_COMPOSITION
        if haskey(sol, k)
            val = sol[k]
            sz = size(val)
            @assert length(sz) == 4 # 3 dims + last index for comps
            ncomp = sz[4]
            x = zeros(ncomp, nc)
            for cell in 1:nc
                I, J, K = cell_ijk(G, cell)
                tot = 0.0
                for i in 1:ncomp
                    v_i = val[I, J, K, i]
                    v_i = max(v_i, ϵ)
                    x[i, cell] = v_i
                    tot += v_i
                end
                err = abs(tot - 1.0)
                @assert err < 0.1 "Too large composition error in cell $cell: sum of compositions $k was $err"
                for i in 1:ncomp
                    x[i, cell] /= tot
                end
            end
        else
            x = missing
        end
        return x
    end

    nph = number_of_phases(sys)
    if is_blackoil
        rs = get_active("RS", to_zero = true)
        rv = get_active("RV", to_zero = true)
    end

    pressure = get_active("PRESSURE")
    s_rem = ones(nc)

    sw = get_active("SWAT")
    if !ismissing(sw)
        s_rem -= sw
    end

    sg = get_active("SGAS")
    if !ismissing(sg)
        s_rem -= sg
    end

    @assert !ismissing(pressure)
    init[:Pressure] = pressure
    if is_blackoil || is_compositional
        if ismissing(sw)
            sw = zeros(nc)
        end
        @assert !ismissing(sg)
        if nph == 3
            init[:ImmiscibleSaturation] = sw
        end
        if is_blackoil
            F_rs = sys.rs_max
            F_rv = sys.rv_max

            so = s_rem
            init[:BlackOilUnknown] = map(
                (w,  o,   g, r,  v, p) -> blackoil_unknown_init(F_rs, F_rv, w, o, g, r, v, p),
                sw, so, sg, rs, rv, pressure)
        else
            @assert !ismissing(sg)
            function get_mole_fraction(k)
                mf = get_active_fraction(k)
                @assert !ismissing(mf)
                return mf
            end
            if haskey(sol, "ZMF")
                z = get_mole_fraction("ZMF")
            else
                x = get_mole_fraction("XMF")
                y = get_mole_fraction("YMF")
                e = 1e-8
                two_ph = count(s -> (s < (1-e) && s > 1e-8), sg)
                if two_ph > 0
                    @warn "XMF/YMF initialization with two-phase conditions does not properly handle multiphase initial conditions. Initial compositions may not be what you expect!"
                end
                z = zeros(size(x))
                for cell in axes(z, 2)
                    for i in axes(z, 1)
                        # TODO: This is wrong for 2ph cells. See warning above.
                        L = sg[i]/(1.0 - sw[i])
                        z[i, cell] = L*x[i] + (1.0-L)*y[i]
                    end
                end
            end
            init[:OverallMoleFractions] = z
        end
    else
        sat = zeros(nph, nc)
        for (idx, phase) in enumerate(get_phases(sys))
            if phase == AqueousPhase()
                sat[idx, :] .= sw
            elseif phase == LiquidPhase()
                sat[idx, :] .= s_rem
            else
                @assert phase == VaporPhase()
                sat[idx, :] .= sg
            end
        end
        init[:Saturations] = sat
    end
    return init
end

function initialize_numerical_aquifers!(init, rmodel, aquifers)
    p = init[:Pressure]
    sys = rmodel.system
    if haskey(init, :OverallMoleFractions)
        # Compositional
        cnames = lowercase.(component_names(sys))
        pos = findfirst(isequal("h2o"), cnames)
        if isnothing(pos)
            pos = findfirst(isequal("water"), cnames)
        end
        if isnothing(pos)
            jutul_message("Did not find water component for aquifers, will not change composition", color = :yellow)
        else
            z = init[:OverallMoleFractions]
            for (id, aquifer) in pairs(aquifers)
                cells = map(x -> x.cell, aquifer.aquifer_cells)
                @. z[:, cells] = 0.0
                @. z[pos, cells] = 1.0
            end
        end
    elseif haskey(init, :ImmiscibleSaturation)
        # Black-oil
        sw = init[:ImmiscibleSaturation]
        for (id, aquifer) in pairs(aquifers)
            cells = map(x -> x.cell, aquifer.aquifer_cells)
            sw[cells] .= 1.0
        end
    else
        # Immiscible model
        ix = findfirst(isequal(AqueousPhase()), get_phases(sys))
        if isnothing(ix)
            s = init[:Saturations]
            for (id, aquifer) in pairs(aquifers)
                cells = map(x -> x.cell, aquifer.aquifer_cells)
                s[:, cells] .= 0.0
                s[ix, cells] .= 1.0
            end
        end
    end

    for (id, aquifer) in pairs(aquifers)
        for aqprm in aquifer.aquifer_cells
            cell = aqprm.cell
            pa = aqprm.pressure
            if isfinite(pa) && pa > DEFAULT_MINIMUM_PRESSURE
                p[cell] = pa
            end
        end
    end
    return init
end

function parse_reservoir(data_file; zcorn_depths = true, repair_zcorn = true, process_pinch = true)
    grid = data_file["GRID"]
    cartdims = grid["cartDims"]
    G = mesh_from_grid_section(grid; repair_zcorn = repair_zcorn, process_pinch = process_pinch)

    # Handle numerical aquifers
    aqunum = get(grid, "AQUNUM", missing)
    aqucon = get(grid, "AQUCON", missing)
    # TODO: Export this properly
    aquifers = GeoEnergyIO.CornerPointGrid.mesh_add_numerical_aquifers!(G, aqunum, aqucon)
    active_ix = G.cell_map
    nc = number_of_cells(G)
    nf = number_of_faces(G)

    permx = grid["PERMX"]
    permy = grid["PERMY"]
    permz = grid["PERMZ"]
    gporo = grid["PORO"]
    gmultpv = get(grid, "MULTPV", Float64[])
    T = promote_type(eltype(permx), eltype(permy), eltype(permz), eltype(gporo), eltype(gmultpv))
    # TODO: PERMYY etc for full tensor perm
    perm = zeros(T, 3, nc)
    poro = zeros(T, nc)
    for (i, c) in enumerate(active_ix)
        perm[1, i] = permx[c]
        perm[2, i] = permy[c]
        perm[3, i] = permz[c]
        poro[i] = gporo[c]
    end
    extra_data_arg = Dict{Symbol, Any}()
    if length(gmultpv) > 0
        multpv = ones(T, nc)
        for (i, c) in enumerate(active_ix)
            v = gmultpv[c]
            if isfinite(v)
                multpv[i] = v
            end
        end
        if !isnothing(aquifers)
            # Avoid MULTPV for aquifer cells
            for (aq_id, aqprm) in pairs(aquifers)
                multpv[aqprm.cell] = 1.0
            end
        end
        extra_data_arg[:pore_volume_multiplier] = multpv
    end




    tranmult = ones(nf)
    tran_override = fill(NaN, nf)
    ijk = map(i -> Jutul.cell_ijk(G, i), 1:nc)
    for (secname, section) in pairs(data_file)
        if haskey(section, "NTG")
            ntg = zeros(nc)
            for (i, c) in enumerate(active_ix)
                ntg[i] = section["NTG"][c]
            end
            extra_data_arg[:net_to_gross] = ntg
        end
        # TODO: This is not 100% robust if edit and grid interact.
        if haskey(section, "PORV")
            extra_data_arg[:pore_volume_override] = section["PORV"][active_ix]
        end
        # Explicit trans given as cell values
        for k in ["TRANX", "TRANY", "TRANZ"]
            if haskey(section, k)
                d = findfirst(isequal(k[end]), ('X', 'Y', 'Z'))
                for (c, val) in enumerate(section[k][active_ix])
                    ijk_c = ijk[c][d]
                    if !isfinite(val) || val < 0.0
                        continue
                    end
                    for face in G.faces.cells_to_faces[c]
                        l, r = G.faces.neighbors[face]
                        if l == c
                            ijk_other = ijk[r][d]
                        else
                            @assert r == c
                            ijk_other = ijk[l][d]
                        end
                        if ijk_other != ijk_c
                            tran_override[face] = val
                        end
                    end
                end
            end
        end

        if haskey(section, "MULTFLT")
            for (fault, vals) in section["MULTFLT"]
                fault_faces = get_mesh_entity_tag(G, Faces(), :faults, Symbol(fault))
                tranmult[fault_faces] *= vals[1]
            end
        end
        mult_keys = ("MULTX", "MULTX-", "MULTY", "MULTY-", "MULTZ", "MULTZ-")
        for k in mult_keys
            if haskey(section, k)
                mult_on_active = section[k][active_ix]
                apply_mult_xyz!(tranmult, k, mult_on_active, G, ijk)
            end
        end
        for k in ["MULTRANX", "MULTRANY", "MULTRANZ"]
            if haskey(section, k)
                fake_mult = ones(cartdims)
                direction = k[end]
                for (I_pair, J_pair, K_pair, val) in section[k]
                    @. fake_mult[
                        I_pair[1]:I_pair[2],
                        J_pair[1]:J_pair[2],
                        K_pair[1]:K_pair[2]
                    ] *= val
                end
                apply_mult_xyz!(tranmult, "MULT$(direction)", fake_mult[active_ix], G, ijk)
            end
        end
    end
    if haskey(grid, "FLUXNUM")
        extra_data_arg[:fluxnum] = grid["FLUXNUM"][active_ix]
    end
    if haskey(grid, "MULTNUM")
        extra_data_arg[:multnum] = grid["MULTNUM"][active_ix]
    end
    if haskey(grid, "OPERNUM")
        extra_data_arg[:opernum] = grid["OPERNUM"][active_ix]
    end
    multregt = get(grid, "MULTREGT", missing)
    if !ismissing(multregt)
        opernum = get(extra_data_arg, :opernum, ones(Int, nc))
        multnum = get(extra_data_arg, :multnum, ones(Int, nc))
        fluxnum = get(extra_data_arg, :fluxnum, ones(Int, nc))
        # Make the table concrete for looping over
        multregt_tab = map(
            x -> (x[1], x[2], x[3], x[4], x[5], only(lowercase(x[6]))),
            multregt
        )
        parser_set_multregt!(tranmult, G, opernum, multnum, fluxnum, multregt_tab, ijk)
    end

    sol = data_file["SOLUTION"]
    has_tempi = haskey(sol, "TEMPI")
    has_tempvd = haskey(sol, "TEMPVD")

    if has_tempvd || has_tempi
        @assert !has_tempvd "TEMPVD not implemented."
        temperature = zeros(nc)
        T_inp = sol["TEMPI"]
        for (i, c) in enumerate(active_ix)
            temperature[i] = convert_to_si(T_inp[i], :Celsius)
        end
        extra_data_arg[:temperature] = temperature
    end
    has(name) = haskey(data_file["RUNSPEC"], name) && data_file["RUNSPEC"][name]

    if has("THERMAL")
        w = has("WATER")
        o = has("OIL")
        g = has("GAS")
        fluid_conductivity = zeros(w+o+g, nc)
        if haskey(grid, "THCONR")
            @assert !haskey(grid, "THCROCK") "THCONR and THCROCK should not be used together."
            extra_data_arg[:rock_thermal_conductivity] = grid["THCONR"][active_ix]
        else
            pos = 1
            for phase in ("WATER", "OIL", "GAS")
                if has(phase)
                    fluid_conductivity[pos, :] .= grid["THC$phase"][active_ix]
                    pos += 1
                end
            end
            if haskey(grid, "THCROCK")
                extra_data_arg[:rock_thermal_conductivity] = grid["THCROCK"][active_ix]
            end
        end
        extra_data_arg[:fluid_thermal_conductivity] = fluid_conductivity
    end

    satnum = GeoEnergyIO.InputParser.get_data_file_cell_region(data_file, :satnum, active = active_ix)
    pvtnum = GeoEnergyIO.InputParser.get_data_file_cell_region(data_file, :pvtnum, active = active_ix)
    eqlnum = GeoEnergyIO.InputParser.get_data_file_cell_region(data_file, :eqlnum, active = active_ix)

    if !isnothing(aquifers)
        for (aq_id, aquifer) in pairs(aquifers)
            for aqprm in aquifer.aquifer_cells
                # Set satnum, pvtnum, static props and verify that cell is present.
                cell = aqprm.cell
                @assert cell <= nc "Numerical aquifer with id $aq_id exceeds number of cells $nc in mesh. Possible failure in aquifer processing."
                satnum[cell] = aqprm.satnum
                pvtnum[cell] = aqprm.pvtnum
                perm[:, cell] .= aqprm.permeability
                poro[cell] = aqprm.porosity
            end
        end
        for (k, v) in extra_data_arg
            if k == :net_to_gross
                for (aq_id, aquifer) in pairs(aquifers)
                    for aqprm in aquifer.aquifer_cells
                        cell = aqprm.cell
                        # Net to gross should not be set for aquifers?
                        v[cell] = 1.0
                    end
                end
            end
        end
    end

    set_scaling_arguments!(extra_data_arg, active_ix, data_file)

    domain = reservoir_domain(G;
        permeability = perm,
        porosity = poro,
        satnum = satnum,
        eqlnum = eqlnum,
        pvtnum = pvtnum,
        pairs(extra_data_arg)...
    )
    if !all(isequal(1.0), tranmult)
        domain[:transmissibility_multiplier, Faces()] = tranmult
    end
    if any(isfinite, tran_override)
        domain[:transmissibility_override, Faces()] = tran_override
    end
    if haskey(grid, "NNC")
        nnc = grid["NNC"]
        if length(nnc) > 0
            domain[:nnc, nothing] = nnc
        end
    end
    edit = get(data_file, "EDIT", Dict())
    if haskey(grid, "DEPTH")
        for (i, c) in enumerate(active_ix)
            domain[:cell_centroids][3, i] = grid["DEPTH"][c]
        end
    elseif haskey(edit, "DEPTH")
        for (i, c) in enumerate(active_ix)
            val = edit["DEPTH"][c]
            if isfinite(val)
                domain[:cell_centroids][3, i] = val
            end
        end
    elseif haskey(grid, "ZCORN") && zcorn_depths
        # Option to use ZCORN points to set depths
        z = get_zcorn_cell_depths(G, grid)
        @. domain[:cell_centroids][3, :] = z
    end

    # pinch_face = get_mesh_entity_tag(G, Faces(), :cpgrid_connection_type, :pinched, throw = false)
    if !isnothing(aquifers)
        domain[:numerical_aquifers, nothing] = aquifers
        vol = domain[:volumes]
        centroids = domain[:cell_centroids]
        for (aq_id, aquifer) in pairs(aquifers)
            for aqprm in aquifer.aquifer_cells
                cell = aqprm.cell
                A = aqprm.area
                L = aqprm.length
                D = aqprm.depth
                vol[cell] = A*L
                centroids[3, cell] = D
            end
        end
    end
    return domain
end

function parser_set_multregt!(tranmult, G, opernum, multnum, fluxnum, multregt, ijk)
    function tsort(x, y)
        if x > y
            return (y, x)
        else
            return (x, y)
        end
    end
    function pair_matchex(pair_kw, pair_reg)
        wildcard1 = pair_kw[1] < 1
        wildcard2 = pair_kw[2] < 1
        if wildcard1 && wilcard2
            return pair_reg[1] != pair_reg[2]
        elseif wildcard1
            return pair_kw[2] == pair_reg[2]
        elseif wildcard2
            return pair_kw[1] == pair_reg[1]
        else
            return pair_kw == pair_reg
        end
    end
    for fno in eachindex(tranmult)
        l, r = G.faces.neighbors[fno]
        opernum_pair = tsort(opernum[l], opernum[r])
        multnum_pair = tsort(multnum[l], multnum[r])
        fluxnum_pair = tsort(fluxnum[l], fluxnum[r])

        for regt in multregt
            pairt = tsort(regt[1], regt[2])
            do_apply = false
            for (pos, coord) in enumerate(('X', 'Y', 'Z'))
                if coord in regt[4] && ijk[l][pos] != ijk[r][pos]
                    do_apply = true
                    break
                end
            end
            if do_apply
                m = regt[3]
                region_type = regt[6]
                if region_type == 'm'
                    pair_to_match = multnum_pair
                elseif region_type == 'o'
                    pair_to_match = opernum_pair
                else
                    @assert region_type == 'f' "Region type was expected to be m, o or f, was $region_type"
                    pair_to_match = fluxnum_pair
                end
                if pair_matchex(pairt, pair_to_match)
                    tranmult[fno] *= m
                end
            end
        end
    end
    return tranmult
end

function apply_mult_xyz!(tranmult, k, mult_on_active, G, ijk)
    nf = length(tranmult)
    if startswith(k, "MULTX")
        pos = 1
    elseif startswith(k, "MULTY")
        pos = 2
    else
        @assert startswith(k, "MULTZ")
        pos = 3
    end
    for fno in 1:nf
        l, r = G.faces.neighbors[fno]
        il = ijk[l][pos]
        ir = ijk[r][pos]
        ok = true
        for i in 1:3
            # Exclude faces that do not match the other logical indices
            # (i.e. NNC or fault related)
            if i == pos
                ok = ok && il != ir
            else
                ok = ok && ijk[l][i] == ijk[r][i]
            end
        end
        if ok
            if il < ir
                low, hi = l, r
            else
                hi, low = l, r
            end
            if endswith(k, "-")
                # Value at current index modifies e.g. K to K-1
                c = hi
            else
                # Value at current index modifies e.g. K to K+1 Pick the
                # lowest cell value (in the sense of the current IJK index)
                # for the face to select a cell-wise MULTX/MULTY/MULTZ for
                # this face.
                c = low
            end
            tranmult[fno] *= mult_on_active[c]
        end
    end
    return tranmult
end

function get_zcorn_cell_depths(g, grid)
    nc = number_of_cells(g)
    z = zeros(nc)
    cartdims = grid["cartDims"]
    zcorn = grid["ZCORN"]
    for c in 1:nc
        i, j, k = cell_ijk(g, c)
        linear_ix = GeoEnergyIO.CornerPointGrid.ijk_to_linear(i, j, k, cartdims)
        get_zcorn(I1, I2, I3) = zcorn[GeoEnergyIO.CornerPointGrid.corner_index(linear_ix, (I1, I2, I3), cartdims)]
        get_pair(I, J) = (get_zcorn(I, J, 0), get_zcorn(I, J, 1))
        l_11, t_11 = get_pair(0, 0)
        l_12, t_12 = get_pair(0, 1)
        l_21, t_21 = get_pair(1, 0)
        l_22, t_22 = get_pair(1, 1)

        z[c] = (1.0/8.0)*(l_11 + t_11 + l_12 + t_12 + l_21 + t_21 + l_22 + t_22)
    end
    return z
end

function set_scaling_arguments!(out, active, data_file)
    function set_scaler!(out, props, phase, prefix)
        phase = uppercase(phase)
        get_scaler_field = x -> get(props, "$(prefix)$x", missing)
        connate = get_scaler_field("S$(phase)L")
        crit = get_scaler_field("S$(phase)CR")
        maxs = get_scaler_field("S$(phase)U")
        if phase == "OW" || phase == "OG"
            maxk = get_scaler_field("KRO")
        else
            maxk = get_scaler_field("KR$phase")
        end
        raw_scalers = (connate, crit, maxs, maxk)
        found = any(!ismissing, raw_scalers)
        nc = length(active)

        if found
            scaler = fill(NaN, (4, nc))
            for (scaler_no, val) in enumerate(raw_scalers)
                if ismissing(val)
                    continue
                end
                for (i, cell) in enumerate(active)
                    scaler[scaler_no, i] = val[cell]
                end
            end
            if prefix == ""
                s = "drainage"
            else
                s = "imbibition"
            end
            out[Symbol("scaler_$(lowercase(phase))_$s")] = scaler
        end
    end

    if haskey(data_file, "PROPS")
        props = data_file["PROPS"]
        for phase in ["W", "OW", "OG", "G"]
            for prefix in ["", "I"]
                set_scaler!(out, props, phase, prefix)
            end
        end
    end
    return out
end

function parse_physics_types(datafile; pvt_region = missing)
    runspec = datafile["RUNSPEC"]
    props = datafile["PROPS"]
    has(name) = haskey(runspec, name) && runspec[name] == true
    has_wat = has("WATER")
    has_oil = has("OIL")
    has_gas = has("GAS")
    has_disgas = has("DISGAS")
    has_vapoil = has("VAPOIL")
    has_thermal = has("THERMAL")
    if has("JUTUL_CO2BRINE")
        # Early termination since this model does not need any PVT specification.
        @assert (has_oil && has_gas) "JUTUL_CO2BRINE currently assumes oil(liquid)-gas systems"
        return (system = :co2brine, pvt = (missing, missing))
    elseif has("CO2STORE")
        error("CO2STORE is not directly supported in JutulDarcy. You can try the converter:\ndata = JutulDarcy.convert_co2store_to_co2_brine(data)")
    end

    is_immiscible = !has_disgas && !has_vapoil
    is_compositional = haskey(runspec, "COMPS")

    phases = []
    rhoS = Vector{Float64}()
    pvt = []
    if haskey(props, "DENSITY")
        deck_densities = props["DENSITY"]
        if deck_densities isa Matrix
            deck_densities = flat_region_expand(deck_density, 3)
        end
        if ismissing(pvt_region)
            if length(deck_densities) > 1
                @warn "Multiple PVT regions found. Picking first one." deck_densities
            end
            pvt_region = 1
        end
        num_pvt_reg = length(deck_densities)
        pvt_region in 1:num_pvt_reg || throw(ArgumentError("Single PVT region found but region $pvt_region was requested."))
        deck_density = deck_densities[pvt_region]

        rhoW_all = map(x -> x[2], deck_densities)
        rhoO_all = map(x -> x[1], deck_densities)
        rhoG_all = map(x -> x[3], deck_densities)

        # Scale densities by dividing by the reference density and multiplying
        # with their PVT region's density to make it consistent when you
        # evaluate density by rhoS*(1/B)
        oil_scale = rhoO_all./rhoO_all[pvt_region]
        gas_scale = rhoG_all./rhoG_all[pvt_region]
        water_scale = rhoW_all./rhoW_all[pvt_region]
        scaling = (
            water_density = water_scale,
            oil_density = oil_scale,
            gas_density = gas_scale,
            rs = gas_scale./oil_scale,
            rv = oil_scale./gas_scale
            )

        rhoOS = deck_density[1]
        rhoWS = deck_density[2]
        rhoGS = deck_density[3]
    else
        @assert is_compositional "DENSITY can only be omitted for compositional models."
        rhoOS = rhoWS = rhoGS = 1.0
        has_oil = true
        has_gas = true
        scaling = missing
    end

    if has_wat
        water_pvt = deck_pvt_water(props, scaling = scaling)
        push!(pvt, water_pvt)
        push!(phases, AqueousPhase())
        push!(rhoS, rhoWS)
    end

    if is_compositional
        push!(pvt, missing)
        push!(phases, LiquidPhase())
        push!(rhoS, rhoOS)
        push!(pvt, missing)
        push!(phases, VaporPhase())
        push!(rhoS, rhoGS)
        if length(props["MW"]) > 1
            jutul_message("EOSNUM", "$(length(props["MW"])) regions active. Only one region supported. Taking the first set of values for all EOS properties.", color = :yellow)
        end

        cnames = copy(props["CNAMES"])
        acf = first(props["ACF"])
        mw = first(props["MW"])
        p_c = first(props["PCRIT"])
        T_c = first(props["TCRIT"])
        if haskey(props, "VCRIT")
            V_c = first(props["VCRIT"])
            if haskey(props, "ZCRIT")
                jutul_message("ZCRIT/VCRIT", "Both VCRIT and ZCRIT found. Using VCRIT for critical volume.", color = :yellow)
            end
        else
            Z_c = first(props["ZCRIT"])
            # J / (mol K)
            # TODO: Double check units
            R = 8.314462618
            V_c = @. Z_c.*R*T_c/p_c
        end

        if haskey(props, "BIC")
            A_ij = first(props["BIC"])
        else
            A_ij = nothing
        end
        mp = MolecularProperty.(mw, p_c, T_c, V_c, acf)
        @assert length(cnames) == length(mp)
        mixture = MultiComponentMixture(mp, A_ij = A_ij, names = cnames)
        if haskey(props, "EOS")
            eos_str = uppercase(first(props["EOS"]))
            if eos_str == "PR"
                eos_type = PengRobinson()
            elseif eos_str == "SRK"
                eos_type = SoaveRedlichKwong()
            elseif eos_str == "RK"
                eos_type = RedlichKwong()
            else
                @assert eos_str == "ZJ" "Unexpected EOS $eos_str: Should be one of PR, SRK, RK, ZJ"
                eos_type = ZudkevitchJoffe()
            end
        else
            eos_type = PengRobinson()
        end
        if haskey(props, "SSHIFT")
            vshift = Tuple(first(props["SSHIFT"]))
        else
            vshift = nothing
        end
        eos = GenericCubicEOS(mixture, eos_type, volume_shift = vshift)
        sys = MultiPhaseCompositionalSystemLV(eos, phases, reference_densities = rhoS)
    else
        if has_oil
            push!(pvt, deck_pvt_oil(props, scaling = scaling))
            push!(phases, LiquidPhase())
            push!(rhoS, rhoOS)
        end

        if has_gas
            push!(pvt, deck_pvt_gas(props, scaling = scaling))
            push!(phases, VaporPhase())
            push!(rhoS, rhoGS)
        end
        sys = pick_system_from_pvt(pvt, rhoS, phases, is_immiscible)
    end
    return (system = sys, pvt = pvt)
end

function pick_system_from_pvt(pvt, rhoS, phases, is_immiscible)
    if is_immiscible
        sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    else
        has_water = length(pvt) == 3
        oil_pvt = pvt[1 + has_water]
        if oil_pvt isa PVTO
            rs_max = map(saturated_table, oil_pvt.tab)
        else
            rs_max = nothing
        end
        gas_pvt = pvt[2 + has_water]
        if gas_pvt isa PVTG
            rv_max = map(saturated_table, gas_pvt.tab)
        else
            rv_max = nothing
        end
        sys = StandardBlackOilSystem(
            rs_max = rs_max,
            rv_max = rv_max,
            phases = phases,
            reference_densities = rhoS
        )
    end
    return sys
end

function data_file_active_tracers(datafile_or_runspec)
    if haskey(datafile_or_runspec, "RUNSPEC")
        rs = datafile_or_runspec["RUNSPEC"]
    else
        rs = datafile_or_runspec
    end
    out = String[]
    possible_tracers = ("POLYMER",)
    for t in possible_tracers
        if haskey(rs, t)
            push!(out, t)
        end
    end
    return out
end

function parse_schedule(domain, sys, datafile; kwarg...)
    schedule = datafile["SCHEDULE"]
    props = datafile["PROPS"]
    return parse_schedule(domain, datafile["RUNSPEC"], props, schedule, sys; kwarg...)
end

function parse_control_steps(runspec, props, schedule, sys)
    rho_s = reference_densities(sys)
    phases = get_phases(sys)

    wells = schedule["WELSPECS"]
    steps = schedule["STEPS"]
    active_tracers = data_file_active_tracers(runspec)
    if length(keys(steps[end])) == 0
        # Prune empty final record
        steps = steps[1:end-1]
    end
    sdict = Dict{String, Any}

    tstep = Vector{Float64}()
    cstep = Vector{Int}()
    well_temp = Dict{String, Float64}()
    compdat = Dict{String, OrderedDict}()
    controls = sdict()
    # "Hidden" well control mirror used with WELOPEN logic
    active_controls = sdict()
    limits = sdict()
    streams = sdict()
    polymer = sdict()
    status = sdict()
    mswell_kw = Dict{String, sdict}()
    function get_and_create_mswell_kw(k::AbstractString, subkey = missing)
        if !haskey(mswell_kw, k)
            mswell_kw[k] = sdict()
        end
        D = mswell_kw[k]
        if ismissing(subkey)
            out = D
        else
            if !haskey(D, subkey)
                D[subkey] = []
            end
            out = D[subkey]
        end
        return out
    end

    well_injection = sdict()
    well_factor = Dict{String, Float64}()
    for k in keys(wells)
        compdat[k] = OrderedDict{NTuple{3, Int}, Any}()
        controls[k] = DisabledControl()
        active_controls[k] = DisabledControl()
        limits[k] = nothing
        streams[k] = nothing
        status[k] = "SHUT"
        polymer[k] = 0.0
        well_injection[k] = nothing
        well_factor[k] = 1.0
        # 0 deg C is strange, but the default for .DATA files.
        well_temp[k] = 273.15 + 0.0
    end
    all_compdat = []
    all_controls = []
    all_limits = []
    all_status = []

    if haskey(runspec, "START")
        start_date = runspec["START"]
    else
        start_date = missing
    end
    current_time = 0.0
    function add_dt!(dt, ctrl_ix)
        if dt ≈ 0
            return
        end
        @assert dt > 0.0 "dt must be positive, attempt to add dt number $(length(tstep)) was $dt at control $ctrl_ix"
        push!(tstep, dt)
        push!(cstep, ctrl_ix)
    end

    skip = ("WELLSTRE", "WINJGAS", "GINJGAS", "GRUPINJE", "WELLINJE", "WEFAC", "WTEMP", "WPIMULT", "WPOLYMER")
    bad_kw = Dict{String, Bool}()
    streams = setup_well_streams()
    for (ctrl_ix, step) in enumerate(steps)
        found_time = false
        streams = parse_well_streams_for_step!(streams, step, props)
        if haskey(step, "WEFAC")
            for wk in step["WEFAC"]
                well_factor[wk[1]] = wk[2]
            end
        end
        if haskey(step, "WTEMP")
            for wk in step["WTEMP"]
                wnm, wt = wk
                well_temp[wnm] = convert_to_si(wt, :Celsius)
            end
        end
        if haskey(step, "WPOLYMER")
            for wk in step["WPOLYMER"]
                name, polymer_concentration = wk
                polymer[name] = polymer_concentration
            end
        end
        for (key, kword) in pairs(step)
            if key == "DATES"
                if ismissing(start_date)
                    @warn "Defaulted date in parsed data and DATES is used. Setting start date to Jan 1. 1983."
                    start_date = DateTime("1983-01-01")
                end
                @assert !found_time
                found_time = true
                for date in kword
                    cdate = start_date + Second(current_time)
                    dt = Float64(Second(date - cdate).value)
                    add_dt!(dt, ctrl_ix)

                    current_time += dt
                end
            elseif key == "TIME"
                @assert !found_time
                found_time = true

                for time in kword
                    dt = time - current_time
                    add_dt!(dt, ctrl_ix)
                    current_time = time
                end
            elseif key == "TSTEP"
                @assert !found_time
                found_time = true
                for dt in kword
                    add_dt!(dt, ctrl_ix)
                    current_time = current_time + dt
                end
            elseif key == "COMPDAT"
                for cd in kword
                    wname, I, J, K1, K2, flag, satnum, WI, diam, Kh, skin, Dfac, dir = cd
                    @assert haskey(wells, wname)
                    head = wells[wname].head
                    if I < 1
                        I = head[1]
                    end
                    if J < 1
                        J = head[2]
                    end
                    entry = compdat[wname]
                    for K in K1:K2
                        perf_mult = 1.0
                        current_completion_index = length(keys(entry)) + 1
                        if haskey(step, "WPIMULT")
                            for wpimult in step["WPIMULT"]
                                wname_wpi, mul_i, I_wpi, J_wpi, K_wpi, start_wpi, end_wpi = wpimult
                                if wname_wpi != wname
                                    continue
                                end
                                if I != I_wpi && I_wpi > 0
                                    continue
                                end
                                if J != J_wpi && J_wpi > 0
                                    continue
                                end
                                if K != K_wpi && K_wpi > 0
                                    continue
                                end
                                if start_wpi > current_completion_index && start_wpi > 0
                                    continue
                                end
                                if end_wpi < current_completion_index && end_wpi > 0
                                    continue
                                end
                                perf_mult *= mul_i
                            end
                        end
                        entry[(I, J, K)] = (
                            open = flag == "OPEN",
                            satnum = satnum,
                            WI = WI,
                            diameter = diam,
                            Kh = Kh,
                            skin = skin,
                            dir = dir,
                            mul = perf_mult,
                            ctrl = ctrl_ix)
                    end
                end
            elseif key in ("WCONINJE", "WCONPROD", "WCONHIST", "WCONINJ", "WCONINJH")
                for wk in kword
                    name = wk[1]
                    controls[name], limits[name], status[name] = keyword_to_control(sys, streams, wk, key,
                        factor = well_factor[name],
                        temperature = well_temp[name],
                        polymer = polymer[name],
                        tracers = active_tracers
                    )
                end
                set_active_controls!(active_controls, controls)
            elseif key == "WELOPEN"
                for wk in kword
                    apply_welopen!(controls, compdat, wk, active_controls)
                end
            elseif key == "WELTARG"
                for wk in kword
                    controls, limits = apply_weltarg!(controls, limits, wk)
                end
                set_active_controls!(active_controls, controls)
            elseif key in skip
                # Already handled
            elseif key in ("WSEGVALV", "COMPSEGS", "WELSEGS")
                for val in kword
                    wname = val[1]
                    ms_storage = get_and_create_mswell_kw(wname, key)
                    for v in val
                        push!(ms_storage, v)
                    end
                end
            else
                bad_kw[key] = true
            end
        end
        if found_time
            push!(all_compdat, deepcopy(compdat))
            push!(all_controls, deepcopy(controls))
            push!(all_limits, deepcopy(limits))
            push!(all_status, deepcopy(status))
        else
            @warn "Did not find supported time kw in step $ctrl_ix: Keys were $(keys(step))."
        end
    end
    for k in keys(bad_kw)
        jutul_message("Unsupported keyword", "Keyword $k was present, but is not supported.", color = :yellow)
    end
    return (
        dt = tstep,
        control_step = cstep,
        controls = all_controls,
        status = all_status,
        completions = all_compdat,
        multisegment = mswell_kw,
        limits = all_limits
    )
end

function set_active_controls!(active_controls, controls)
    for (name, ctrl) in pairs(controls)
        if ctrl isa DisabledControl
            continue
        end
        active_controls[name] = ctrl
    end
    return active_controls
end

function setup_well_streams()
    return (streams = Dict{String, Any}(), wells = Dict{String, String}())
end

function parse_well_streams_for_step!(streams, step, props)
    well_streams = streams.wells
    streams = streams.streams
    if haskey(props, "STCOND")
        std = props["STCOND"]
        T = convert_to_si(std[1], :Celsius)
        p = std[2]
    else
        @debug "Defaulted STCOND..."
        T = convert_to_si(15.56, :Celsius)
        p = convert_to_si(1.0, :atm)
    end

    if haskey(step, "WELLSTRE")
        for stream in step["WELLSTRE"]
            mix = Float64.(stream[2:end])
            streams[first(stream)] = (mole_fractions = mix, cond = (p = p, T = T))
        end
    end
    if haskey(step, "WINJGAS")
        winjgas = step["WINJGAS"]
        for winjgas in step["WINJGAS"]
            @assert uppercase(winjgas[2]) == "STREAM"
            well_streams[winjgas[1]] = winjgas[3]
        end
    end
    for kw in ["GINJGAS", "GRUPINJE", "WELLINJE"]
        if haskey(step, kw)
            @warn "Stream keyword $kw was found but is not supported."
        end
    end

    return (streams = streams, wells = well_streams)
end

function keyword_to_control(sys, streams, kw, k::String; kwarg...)
    return keyword_to_control(sys, streams, kw, Val(Symbol(k)); kwarg...)
end

function keyword_to_control(sys, streams, kw, ::Val{:WCONPROD}; kwarg...)
    rho_s = reference_densities(sys)
    phases = get_phases(sys)

    flag = kw[2]
    ctrl = kw[3]
    orat = kw[4]
    wrat = kw[5]
    grat = kw[6]
    lrat = kw[7]
    resv = kw[8]
    bhp = kw[9]
    return producer_control(sys, flag, ctrl, orat, wrat, grat, lrat, bhp; resv = resv, kwarg...)
end

function keyword_to_control(sys, streams, kw, ::Val{:WCONHIST}; kwarg...)
    rho_s = reference_densities(sys)
    phases = get_phases(sys)
    # 1 name
    flag = kw[2]
    ctrl = kw[3]
    orat = kw[4]
    wrat = kw[5]
    grat = kw[6]
    lrat = wrat + orat
    # TODO: Pass thp on
    thp = kw[9]
    bhp = kw[10]

    return producer_control(sys, flag, ctrl, orat, wrat, grat, lrat, bhp; is_hist = true, kwarg...)
end

function producer_limits(; bhp = Inf, lrat = Inf, orat = Inf, wrat = Inf, grat = Inf)
    lims = Dict{Symbol, Any}()
    if isfinite(bhp)
        lims[:bhp] = bhp
    end
    if isfinite(lrat)
        lims[:lrat] = -lrat
    end
    if isfinite(orat)
        lims[:orat] = -orat
    end
    if isfinite(wrat)
        lims[:wrat] = -wrat
    end
    if isfinite(grat)
        lims[:grat] = -grat
    end
    return NamedTuple(pairs(lims))
end

function producer_control(sys, flag, ctrl, orat, wrat, grat, lrat, bhp; is_hist = false, temperature = NaN, resv = 0.0, tracers = String[], polymer = 0.0, kwarg...)
    rho_s = reference_densities(sys)
    phases = get_phases(sys)

    if flag == "SHUT" || flag == "STOP"
        ctrl = DisabledControl()
        lims = nothing
    else
        is_rate = true
        @assert flag == "OPEN"
        if ctrl == "LRAT"
            self_val = -lrat
            t = SurfaceLiquidRateTarget(self_val)
        elseif ctrl == "WRAT"
            self_val = -wrat
            t = SurfaceWaterRateTarget(self_val)
        elseif ctrl == "ORAT"
            self_val = -orat
            t = SurfaceOilRateTarget(self_val)
        elseif ctrl == "GRAT"
            self_val = -grat
            t = SurfaceGasRateTarget(self_val)
        elseif ctrl == "BHP"
            self_val = bhp
            t = BottomHolePressureTarget(self_val)
            is_rate = false
        elseif ctrl == "RESV"
            if is_hist
                self_val = -(wrat + orat + grat)
                w = [wrat, orat, grat]
                w = w./sum(w)
                t = HistoricalReservoirVoidageTarget(self_val, w)
            else
                self_val = resv
                w = [1.0, 1.0, 1.0]
                t = ReservoirVoidageTarget(self_val, w)
            end
        else
            error("$ctype control not supported")
        end
        if is_rate && abs(self_val) < MIN_ACTIVE_WELL_RATE
            @debug "Producer with $ctrl disabled due to zero rate." abs(self_val)
            ctrl = DisabledControl()
        else
            ctrl = ProducerControl(t; kwarg...)
        end
        if is_hist
            self_symbol = translate_target_to_symbol(t, shortname = true)
            # Put pressure slightly above 1 atm to avoid hard limit.
            if self_symbol == :resv_history
                lims = (; :bhp => 1.01*si_unit(:atm))
            else
                lims = (; :bhp => 1.01*si_unit(:atm), self_symbol => self_val)
            end
        else
            lims = producer_limits(bhp = bhp, orat = orat, wrat = wrat, grat = grat, lrat = lrat)
        end
    end
    return (ctrl, lims, flag)
end

function injector_limits(; bhp = Inf, surface_rate = Inf, reservoir_rate = Inf)
    lims = Dict{Symbol, Any}()
    if isfinite(bhp)
        lims[:bhp] = bhp
    end
    if isfinite(surface_rate)
        lims[:rate] = surface_rate
    end
    if isfinite(reservoir_rate)
        lims[:resv_rate] = reservoir_rate
    end
    return NamedTuple(pairs(lims))
end

function injector_control(sys, streams, name, flag, type, ctype, surf_rate, res_rate, bhp; is_hist = false, tracers = String[], polymer = 0.0, kwarg...)
    if occursin('*', flag)
        # This is a bit of a hack.
        flag = "OPEN"
    end
    well_is_disabled = flag == "SHUT" || flag == "STOP"
    if !well_is_disabled
        @assert flag == "OPEN" "Unsupported well flag: $flag"
        if ctype == "RATE"
            is_rate = true
            if isfinite(surf_rate)
                t = TotalRateTarget(surf_rate)
            else
                well_is_disabled = true
            end
        elseif ctype == "BHP"
            is_rate = false
            if isfinite(bhp)
                t = BottomHolePressureTarget(bhp)
            else
                well_is_disabled = true
            end
        elseif ctype == "RESV"
            is_rate = false
            if isfinite(res_rate) && res_rate > MIN_ACTIVE_WELL_RATE
                t = TotalReservoirRateTarget(res_rate)
            else
                well_is_disabled = true
            end
        else
            # RESV, GRUP, THP
            error("$ctype control not supported")
        end
        rho, mix, phases_mix = select_injector_mixture_spec(sys, name, streams, type)
        if is_rate && surf_rate < MIN_ACTIVE_WELL_RATE
            @debug "Disabling injector $name with $ctype ctrl due to zero rate" surf_rate
            well_is_disabled = true
        end
    end
    if well_is_disabled
        ctrl = DisabledControl()
        lims = nothing
    else
        ntracers = length(tracers)
        if ntracers > 0
            tracer_vals = zeros(ntracers)
            for (i, t) in enumerate(tracers)
                if t == "POLYMER"
                    tracer_vals[i] = polymer
                else
                    error("Tracer $t not supported.")
                end
            end
        else
            tracer_vals = missing
        end
        ctrl = InjectorControl(t, mix; density = rho, phases = phases_mix, tracers = tracer_vals, kwarg...)
        if is_hist
            # TODO: This magic number comes from MRST.
            bhp_lim = convert_to_si(6895.0, :bar)
        else
            bhp_lim = bhp
        end
        lims = injector_limits(bhp = bhp_lim, surface_rate = surf_rate, reservoir_rate = res_rate)
    end
    return (ctrl, lims, flag)
end

function select_injector_mixture_spec(sys::Union{ImmiscibleSystem, StandardBlackOilSystem}, name, streams, type)
    rho_s = reference_densities(sys)
    phases = get_phases(sys)
    mix = Float64[]
    rho = 0.0
    nph = number_of_phases(sys)
    phases_mix = Tuple{Int, Float64}[]
    ix = 1
    for (phase, rho_ph) in zip(phases, rho_s)
        if phase == LiquidPhase()
            v = Float64(type == "OIL")
        elseif phase == AqueousPhase()
            v = Float64(type == "WATER" || type == "WAT")
        else
            @assert phase isa VaporPhase
            v = Float64(type == "GAS")
        end
        push!(mix, v)
        push!(phases_mix, (ix, v))
        rho += rho_ph*v
        ix += 1
    end
    @assert sum(mix) ≈ 1.0 "Expected mixture to sum to 1, was $mix for type $type (declared phases: $phases)"
    return (rho, mix, phases_mix)
end

function select_injector_mixture_spec(sys::CompositionalSystem, name, streams, type)
    eos = sys.equation_of_state
    props = eos.mixture.properties
    rho_s = reference_densities(sys)
    phases = get_phases(sys)
    has_water = has_other_phase(sys)
    offset = Int(has_water)
    ncomp = number_of_components(sys)
    # Well stream will be molar fractions.
    mix = zeros(Float64, ncomp)
    if uppercase(type) == "WATER" || uppercase(type) == "WAT"
        has_water || throw(ArgumentError("Cannot have WATER injector without water phase."))
        mix[end] = 1.0
        rho = rho_s[1]
        phases_mix = ((1, 1.0), (2, 0.0), (3, 0.0))
    else
        if !haskey(streams.wells, name)
            throw(ArgumentError("Well $name does not have a stream declared with well type $type."))
        end
        stream_id = streams.wells[name]
        stream = streams.streams[stream_id]

        ϵ = MultiComponentFlash.MINIMUM_COMPOSITION
        z = map(z_i -> max(z_i, ϵ), stream.mole_fractions)
        z /= sum(z)
        cond = stream.cond

        z_mass = map(
            (z_i, prop) -> max(z_i*prop.mw, ϵ),
            z, props
        )
        z_mass /= sum(z_mass)
        for i in eachindex(z_mass)
            mix[i] = z_mass[i]
        end
        @assert sum(mix) ≈ 1.0 "Sum of mixture was $(sum(mix)) != 1 for mole mixture $(z) as mass $z_mass"

        flash_cond = (p = cond.p, T = cond.T, z = z)
        flash = MultiComponentFlash.flashed_mixture_2ph(eos, flash_cond)
        rho_l, rho_v = MultiComponentFlash.mass_densities(eos, cond.p, cond.T, flash)
        S_l, S_v = MultiComponentFlash.phase_saturations(eos, cond.p, cond.T, flash)
        rho = S_l*rho_l + S_v*rho_v
        if uppercase(type) == "GAS"
            sg = 1.0
        else
            sg = 0.0
        end
        if has_water
            phases_mix = ((1, 1.0), (2, 1.0-sg), (3, sg))
        else
            phases_mix = ((1, 1.0-sg), (2, sg))
        end
    end
    return (rho, mix, phases_mix)
end

function select_injector_mixture_spec(sys::Symbol, name, streams, type)
    return select_injector_mixture_spec(Val(sys), name, streams, type)
end

function keyword_to_control(sys, streams, kw, ::Val{:WCONINJE}; kwarg...)
    # TODO: Expand to handle mixture etc.
    name = kw[1]
    type = kw[2]
    flag = kw[3]
    ctype = kw[4]
    surf_rate = kw[5]
    res_rate = kw[6]
    bhp = kw[7]
    return injector_control(sys, streams, name, flag, type, ctype, surf_rate, res_rate, bhp; kwarg...)
end

function keyword_to_control(sys, streams, kw, ::Val{:WCONINJH}; kwarg...)
    name = kw[1]
    type = kw[2]
    flag = kw[3]
    surf_rate = kw[4]
    bhp = kw[5]
    ctype = kw[12]
    # TODO: Expand to handle mixture etc.
    res_rate = Inf
    return injector_control(sys, streams, name, flag, type, ctype, surf_rate, res_rate, bhp, is_hist = true; kwarg...)
end

function well_completion_sortperm(domain, wspec, order_t0, wc, dir)
    order_t = lowercase(order_t0)
    @assert order_t in ("track", "input", "depth") "Invalid order for well: $order_t0"
    centroid(dim) = domain[:cell_centroids][dim, wc]
    n = length(wc)
    if n < 2 || order_t == "input"
        # Do nothing.
        sorted = eachindex(wc)
    elseif order_t == "depth" || all(isequal("Z"), dir)
        z = centroid(3)
        sorted = sortperm(z)
    else
        sorted = Int[]
        @assert order_t == "track"
        x = centroid(1)
        y = centroid(2)
        z = centroid(3)
        # Make copies so we can safely remove values as we go.
        wc = copy(wc)
        original_ix = collect(1:n)
        dir = lowercase.(copy(dir))
        g = physical_representation(domain)
        ijk = map(ix -> cell_ijk(g, ix), wc)

        function remove_candidate!(ix)
            deleteat!(x, ix)
            deleteat!(y, ix)
            deleteat!(z, ix)
            deleteat!(wc, ix)
            deleteat!(dir, ix)
            deleteat!(ijk, ix)
            deleteat!(original_ix, ix)
        end
        function add_to_sorted!(ix)
            @assert ix > 0 && ix <= length(wc) "Algorithm failure. Programming error?"
            push!(sorted, original_ix[ix])
            previous_ix = closest_ix
            previous_coord = (x[ix], y[ix], z[ix])
            previous_ijk = ijk[ix]
            previous_dir = dir[ix]
            remove_candidate!(ix)
            return (previous_ix, previous_coord, previous_ijk, previous_dir)
        end

        # Pick closest cell to head to start with
        closest_ix = 0
        closest_ij_distance = typemax(Int)
        lowest_z = Inf
        I_head, J_head = wspec.head
        for (i, c) in enumerate(wc)
            z_i = z[i]
            I, J, = ijk[i]
            d = abs(I - I_head) + abs(J - J_head)
            if d == closest_ij_distance
                new_minimum = z_i < lowest_z
            elseif d < closest_ij_distance
                new_minimum = true
            else
                new_minimum = false
            end

            if new_minimum
                closest_ij_distance = d
                closest_ix = i
                lowest_z = z_i
            end
        end
        prev_ix, prev_coord, prev_ijk, prev_dir = add_to_sorted!(closest_ix)
        start = wspec.head
        use_dir = true
        while length(wc) > 0
            closest_ix = 0
            closest_xyz_distance = Inf
            if use_dir
                closest_ijk_distance = typemax(Int)
                if prev_dir == "x"
                    dim = 1
                    coord = x
                elseif prev_dir == "y"
                    dim = 2
                    coord = y
                else
                    @assert prev_dir == "z"
                    dim = 3
                    coord = z
                end
            else
                closest_ijk_distance = Inf
            end
            for (i, c) in enumerate(wc)
                if use_dir
                    d_ijk = abs(ijk[i][dim] - prev_ijk[dim])
                    d_xyz = abs(coord[i] - prev_coord[dim])
                    if d_ijk == closest_ijk_distance
                        new_minimum = d_xyz < closest_xyz_distance
                    elseif d_ijk < closest_ijk_distance
                        new_minimum = true
                    else
                        new_minimum = false
                    end
                else
                    coord = (x[i], y[i], z[i])
                    d_ijk = norm(ijk[i] .- prev_ijk, 2)
                    d_xyz = norm(coord .- prev_coord, 2)
                    # d_xyz = abs(coord[i] - prev_coord[dim])
                    if d_xyz == closest_xyz_distance
                        new_minimum = d_ijk < closest_ijk_distance
                    elseif d_xyz < closest_xyz_distance
                        new_minimum = true
                    else
                        new_minimum = false
                    end
                end
                if new_minimum
                    closest_ix = i
                    closest_ijk_distance = d_ijk
                    closest_xyz_distance = d_xyz
                end
            end
            prev_ix, prev_coord, prev_ijk, prev_dir = add_to_sorted!(closest_ix)
        end
    end
    @assert sort(sorted) == 1:n "$sorted was not $(1:n)"
    @assert length(sorted) == n
    return sorted
end

function apply_weltarg!(controls, limits, wk)
    well, ctype, value = wk
    ctype = lowercase(ctype)
    ctrl = controls[well]
    if ctrl isa DisabledControl
        jutul_message("$well control", "Cannot use WELTARG $ctype = $value on disabled well $well.", color = :red)
    else
        is_injector = ctrl isa InjectorControl
        limit = limits[well]
        limit = OrderedDict{Symbol, Float64}(pairs(limit))
        if ctype == "bhp"
            new_target = BottomHolePressureTarget(value)
            limit[Symbol(ctype)] = value
        else
            if is_injector
                rate_target = value
            else
                rate_target = -value
            end
            if ctype == "orat"
                new_target = SurfaceOilRateTarget(rate_target)
            elseif ctype == "grat"
                new_target = SurfaceGasRateTarget(rate_target)
            elseif ctype == "wrat"
                new_target = SurfaceWaterRateTarget(rate_target)
            elseif ctype == "rate"
                new_target = TotalRateTarget(rate_target)
            elseif ctype == "lrat"
                new_target = SurfaceLiquidRateTarget(rate_target)
            elseif ctype == "resv"
                new_target = TotalReservoirRateTarget(rate_target)
            else
                error("WELTARG $ctype is not yet supported.")
            end
            limit[Symbol(ctype)] = rate_target
        end
        limits[well] = convert_to_immutable_storage(limit)
        controls[well] = replace_target(ctrl, new_target)
    end
    return (controls, limits)
end

function apply_welopen!(controls, compdat, wk, controls_if_active)
    name, flag, I, J, K, first_num, last_num = wk
    # TODO: Handle shut in a better way
    flag = uppercase(flag)
    @assert flag in ("OPEN", "SHUT", "STOP")
    is_open = flag == "OPEN"
    if I == J == K == first_num == last_num == -1
        # Applies to well
        if is_open
            controls[name] = controls_if_active[name]
        else
            controls[name] = DisabledControl()
        end
    else
        cdat = compdat[name]
        ijk = collect(keys(cdat))
        nperf = length(ijk)
        first_num = max(first_num, 1)
        if last_num < 1
            last_num = nperf
        end
        for i in 1:nperf
            if i < first_num
                continue
            elseif i > last_num
                continue
            else
                I_i, J_i, K_i = ijk[i]
                current_cdat = cdat[ijk[i]]
                is_match(ix, ix_i) = ix < 1 || ix_i == ix
                if is_match(I, I_i) && is_match(J, J_i) && is_match(K, K_i)
                    c = OrderedDict{Symbol, Any}()
                    for (k, v) in pairs(current_cdat)
                        c[k] = v
                    end
                    c[:open] = is_open
                    cdat[ijk[i]] = (; c...)
                end
            end
        end
    end
end

function parse_aquifer_bc(model, datafile, sys)
    bc = nothing
    sources = nothing
    sol = datafile["SOLUTION"]
    aquchwat = get(sol, "AQUCHWAT", nothing)
    reservoir = reservoir_domain(model)
    g = physical_representation(reservoir)
    grav = gravity_constant
    if isnothing(aquchwat)
        bc = nothing
    else
        rho, mix, = select_injector_mixture_spec(sys, "AQUIFER", missing, "WATER")
        mix = tuple(mix...)
        bc = FlowBoundaryCondition{Int, Float64, typeof(mix)}[]
        acon = sol["AQUANCON"]
        aquifers = unique(map(first, acon))
        for aqno in aquifers
            aquifer = aquchwat[aqno]
            datum = aquifer[2]
            bc_type = aquifer[3]
            datum_val = aquifer[4]
            if bc_type == "HEAD"
                p_datum = si_unit(:atmosphere) + grav*rho*(datum_val - datum)
            else
                p_datum = datum_val
            end
            PI = aquifer[5]
            minv, maxv = aquifer[8:9]
            if isfinite(minv) || isfinite(maxv)
                jutul_message("Aquifer $aqno", "Non-defaulted pressure limits ($minv,$maxv) detected. They will be ignored.")
            end
            ignore_depth = aquifer[13] == "YES"
            temperature = aquifer[15]
            if isnan(temperature)
                temperature = 273.15
            end
            # Areas, cells, faces of candidates
            areas = Float64[]
            multipliers = Float64[]
            cells = Int[]
            for con in acon
                if con[1] != aqno
                    continue
                end
                multiplier = con[10]
                area = con[9]
                bfaces, bcells = aquifer_connections_to_boundary_faces(g, con[2], con[3], con[4], con[5], con[6], con[7], con[8])
                for c in bcells
                    push!(cells, c)
                    push!(multipliers, multiplier)
                end
                if isnan(area)
                    for f in bfaces
                        push!(areas, reservoir[:boundary_areas][f])
                    end
                else
                    for _ in bfaces
                        push!(areas, area)
                    end
                end
            end
            area_total = sum(areas)
            @. areas /= area_total
            # Scale by typical water viscosity
            scale = 0.01*si_unit(:poise)
            PI_actual = PI.*scale.*areas.*multipliers
            cc = reservoir[:cell_centroids]
            for (c, PI_c) in zip(cells, PI_actual)
                if ignore_depth
                    p = p_datum
                else
                    depth = cc[3, c]
                    p = p_datum + grav*rho*(depth - datum)
                end
                bc_c = FlowBoundaryCondition(c, p, temperature,
                    trans_flow = PI_c,
                    trans_thermal = 0.0, # Maybe not correct?
                    density = rho,
                    fractional_flow = mix
                )
                push!(bc, bc_c)
            end
        end
        if length(bc) == 0
            bc = nothing
        end
    end
    # TODO: Check that there are no aquifers in schedule since this is not
    # handled.

    return (bc, sources)
end

function aquifer_connections_to_boundary_faces(g, Imin, Imax, Jmin, Jmax, Kmin, Kmax, type)
    if Imax == -1
        Imax = typemax(Int)
    end
    if Jmax == -1
        Jmax = typemax(Int)
    end
    if Kmax == -1
        Kmax = typemax(Int)
    end
    if type in ("I+", "X+")
        t = :right
    elseif type in ("I-", "X-")
        t = :left
    elseif type in ("J+", "Y+")
        t = :upper
    elseif type in ("J-", "Y-")
        t = :lower
    elseif type in ("K+", "Z+")
        t = :bottom
    elseif type in ("K-", "Z-")
        t = :top
    else
        error("Bad face type $type")
    end
    candidates = get_mesh_entity_tag(g, BoundaryFaces(), :direction, t)
    cells = g.boundary_faces.neighbors[candidates]
    bfaces = Int[]
    bcells = Int[]

    for (f, c) in zip(candidates, cells)
        i, j, k = cell_ijk(g, c)
        if i < Imin || i > Imax
            continue
        end
        if j < Jmin || j > Jmax
            continue
        end
        if k < Kmin || k > Kmax
            continue
        end
        push!(bfaces, f)
        push!(bcells, c)
    end
    return (bfaces, bcells)
end
