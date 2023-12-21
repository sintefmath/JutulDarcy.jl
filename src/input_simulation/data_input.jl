export simulate_data_file, setup_case_from_data_file

"""
    simulate_data_file(inp; parse_arg = NamedTuple(), kwarg...)

Simulate standard input file (with extension .DATA). `inp` can either be the
output from `setup_case_from_parsed_data` or a String for the path of an input file.

Additional arguments are passed onto `simulate_reservoir`. Extra inputs to the
parser can be sent as a `parse_arg` `NamedTuple`.
"""
function simulate_data_file(data; setup_arg = NamedTuple(), kwarg...)
    case = setup_case_from_data_file(data; setup_arg...)
    result = simulate_reservoir(case; kwarg...)
    result.extra[:case] = case
    return result
end

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

function setup_case_from_parsed_data(datafile; simple_well = true, kwarg...)
    sys, pvt = parse_physics_types(datafile)
    is_blackoil = sys isa StandardBlackOilSystem
    is_compositional = sys isa CompositionalSystem
    domain = parse_reservoir(datafile)
    wells, controls, limits, cstep, dt, well_mul = parse_schedule(domain, sys, datafile; simple_well = simple_well)

    model, parameters0 = setup_reservoir_model(domain, sys; wells = wells, kwarg...)
    for (k, submodel) in pairs(model.models)
        if submodel.system isa MultiPhaseSystem
            # Modify secondary variables
            if !is_compositional
                svar = submodel.secondary_variables
                # PVT
                pvt = tuple(pvt...)
                rho = DeckDensity(pvt)
                if sys isa StandardBlackOilSystem
                    set_secondary_variables!(submodel, ShrinkageFactors = JutulDarcy.DeckShrinkageFactors(pvt))
                end
                mu = DeckViscosity(pvt)
                set_secondary_variables!(submodel, PhaseViscosities = mu, PhaseMassDensities = rho)
            end
            if k == :Reservoir
                rs = datafile["RUNSPEC"]
                oil = haskey(rs, "OIL")
                water = haskey(rs, "WATER")
                gas = haskey(rs, "GAS")

                JutulDarcy.set_deck_specialization!(submodel, datafile["PROPS"], domain[:satnum], oil, water, gas)
            end
        end
    end
    parameters = setup_parameters(model)
    if haskey(datafile["PROPS"], "SWL")
        G = physical_representation(domain)
        swl = vec(datafile["PROPS"]["SWL"])
        parameters[:Reservoir][:ConnateWter] .= swl[G.cell_map]
    end
    forces = parse_forces(model, wells, controls, limits, cstep, dt, well_mul)
    state0 = parse_state0(model, datafile)
    return JutulCase(model, dt, forces, state0 = state0, parameters = parameters)
end

function parse_well_from_compdat(domain, wname, v, wspecs; simple_well = true)
    wc, WI, open = compdat_to_connection_factors(domain, v)
    rd = wspecs.ref_depth
    if isnan(rd)
        rd = nothing
    end
    W = setup_well(domain, wc, name = Symbol(wname), WI = WI, reference_depth = rd, simple_well = simple_well)
    return (W, wc, WI, open)
end

function compdat_to_connection_factors(domain, v)
    G = physical_representation(domain)
    K = domain[:permeability]

    ij_ix = collect(keys(v))
    wc = map(i -> cell_index(G, i), ij_ix)

    getf(k) = map(x -> v[x][k], ij_ix)
    open = getf(:open)
    d = getf(:diameter)
    Kh = getf(:Kh)
    WI = getf(:WI)
    skin = getf(:skin)
    dir = getf(:dir)

    for i in eachindex(WI)
        W_i = WI[i]
        if isnan(W_i) || W_i < 0.0
            c = wc[i]
            if K isa AbstractVector
                k_i = K[c]
            else
                k_i = K[:, c]
            end
            WI[i] = compute_peaceman_index(G, k_i, d[i]/2, c, skin = skin[i], Kh = Kh[i], dir = Symbol(dir[i]))
        end
    end
    return (wc, WI, open)
end

function parse_schedule(domain, runspec, props, schedule, sys; simple_well = true)
    G = physical_representation(domain)

    dt, cstep, controls, completions, limits = parse_control_steps(runspec, props, schedule, sys)
    completions, bad_wells = filter_inactive_completions!(completions, G)
    handle_wells_without_active_perforations!(bad_wells, completions, controls, limits)
    ncomp = length(completions)
    wells = []
    well_mul = []
    for (k, v) in pairs(completions[end])
        W, wc_base, WI_base, open = parse_well_from_compdat(domain, k, v, schedule["WELSPECS"][k]; simple_well = simple_well)
        wi_mul = zeros(length(WI_base), ncomp)
        for (i, c) in enumerate(completions)
            compdat = c[k]
            wc, WI, open = compdat_to_connection_factors(domain, compdat)
            for (c, wi, is_open) in zip(wc, WI, open)
                compl_idx = findfirst(isequal(c), wc_base)
                if isnothing(compl_idx)
                    # This perforation is missing, leave at zero.
                    continue
                end
                wi_mul[compl_idx, i] = wi*is_open/WI_base[compl_idx]
            end
            # @. wi_mul[:, i] = (WI*open)/WI_base
        end
        push!(wells, W)
        if all(isapprox(1.0), wi_mul)
            push!(well_mul, nothing)
        else
            push!(well_mul, wi_mul)
        end
    end
    return (wells, controls, limits, cstep, dt, well_mul)
end

function filter_inactive_completions!(completions_vector, g)
    tuple_active = Dict{Tuple{String, NTuple{3, Int}}, Bool}()
    for completions in completions_vector
        for (well, completion_set) in pairs(completions)
            for tupl in keys(completion_set)
                k = (well, tupl)
                if !haskey(tuple_active, k)
                    ix = cell_index(g, tupl, throw = false)
                    if isnothing(ix)
                        @warn "Well $well completion $tupl was declared using COMPDAT, but that cell is not active in model. Skipped."
                        tuple_active[k] = false
                    else
                        tuple_active[k] = true
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

function parse_forces(model, wells, controls, limits, cstep, dt, well_mul)
    forces = []
    for (ctrl, lim) in zip(controls, limits)
        ctrl_s = Dict{Symbol, Any}()
        for (k, v) in pairs(ctrl)
            ctrl_s[Symbol(k)] = v
        end
        lim_s = Dict{Symbol, Any}()
        for (k, v) in pairs(lim)
            lim_s[Symbol(k)] = v
        end
        f = setup_reservoir_forces(model, control = ctrl_s, limits = lim_s)
        push!(forces, f)
    end
    return forces[cstep]
end

function parse_state0(model, datafile)
    rmodel = JutulDarcy.reservoir_model(model)
    init = Dict{Symbol, Any}()
    sol = datafile["SOLUTION"]

    if haskey(sol, "EQUIL")
        init = parse_state0_equil(rmodel, datafile)
    else
        init = parse_state0_direct_assignment(rmodel, datafile)
    end
    state0 = setup_reservoir_state(model, init)
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
                (w,  o,   g, r,  v, p) -> JutulDarcy.blackoil_unknown_init(F_rs, F_rv, w, o, g, r, v, p),
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

function parse_reservoir(data_file)
    grid = data_file["GRID"]
    cartdims = grid["cartDims"]
    actnum = get_effective_actnum(grid)
    if haskey(grid, "COORD")
        coord = grid["COORD"]
        zcorn = grid["ZCORN"]
        primitives = JutulDarcy.cpgrid_primitives(coord, zcorn, cartdims, actnum = actnum)
        G = JutulDarcy.grid_from_primitives(primitives)
    else
        @assert haskey(grid, "DX")
        @assert haskey(grid, "DY")
        @assert haskey(grid, "DZ")
        @assert haskey(grid, "TOPS")
        @warn "DX+DY+DZ+TOPS format is only supported if all cells are equally sized. If you get an error, this is the cause."
        @assert all(actnum)
        dx = only(unique(grid["DX"]))
        dy = only(unique(grid["DY"]))
        dz = only(unique(grid["DZ"]))
        tops = only(unique(grid["TOPS"]))
        G = CartesianMesh(cartdims, cartdims.*(dx, dy, dz))
        # We always want to return an unstructured mesh.
        G = UnstructuredMesh(G)
    end
    active_ix = G.cell_map
    nc = number_of_cells(G)
    # TODO: PERMYY etc for full tensor perm
    perm = zeros(3, nc)
    poro = zeros(nc)
    for (i, c) in enumerate(active_ix)
        perm[1, i] = grid["PERMX"][c]
        perm[2, i] = grid["PERMY"][c]
        perm[3, i] = grid["PERMZ"][c]
        poro[i] = grid["PORO"][c]
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
        temp_arg = (temperature = temperature, )
    else
        temp_arg = NamedTuple()
    end
    satnum = JutulDarcy.InputParser.table_region(data_file, :saturation, active = active_ix)
    eqlnum = JutulDarcy.InputParser.table_region(data_file, :equil, active = active_ix)

    domain = reservoir_domain(G;
        permeability = perm,
        porosity = poro,
        satnum = satnum,
        eqlnum = eqlnum,
        temp_arg...
    )
end

function get_effective_actnum(g)
    if haskey(g, "ACTNUM")
        actnum = copy(g["ACTNUM"])
    else
        actnum = fill(true, prod(grid["cartDims"]))
    end
    handle_zero_effective_porosity!(actnum, g)
    return actnum
end

function handle_zero_effective_porosity!(actnum, g)
    if haskey(g, "MINPV")
        minpv = g["MINPV"]
    else
        minpv = 1e-6
    end
    added = 0
    active = 0

    if haskey(g, "PORV")
        porv = G["PORV"]
        for i in eachindex(actnum)
            if actnum[i]
                pv = porv[i]
                active += active
                if pv < minpv
                    added += 1
                    actnum[i] = false
                end
            end
        end
    elseif haskey(g, "PORO")
        if haskey(g, "ZCORN")
            zcorn = g["ZCORN"]
            coord = reshape(g["COORD"], 6, :)'
            cartdims = g["cartDims"]
        else
            zcorn = coord = cartdims = missing
        end
        # Have to handle zero or negligble porosity.
        if haskey(g, "NTG")
            ntg = g["NTG"]
        else
            ntg = ones(size(actnum))
        end
        poro = g["PORO"]
        for i in eachindex(actnum)
            if actnum[i]
                vol = zcorn_volume(g, zcorn, coord, cartdims, i)
                pv = poro[i]*ntg[i]*vol
                active += active
                if pv < minpv
                    added += 1
                    actnum[i] = false
                end
            end
        end
    end
    @debug "$added disabled cells out of $(length(actnum)) due to low effective pore-volume."
    return actnum
end

function zcorn_volume(g, zcorn, coord, dims, linear_ix)
    if ismissing(zcorn)
        return 1.0
    end
    nx, ny, nz = dims
    i, j, k = linear_to_ijk(linear_ix, dims)

    get_zcorn(I1, I2, I3) = zcorn[corner_index(linear_ix, (I1, I2, I3), dims)]
    get_pair(I, J) = (get_zcorn(I, J, 0), get_zcorn(I, J, 1))
    function pillar_line(I, J)
        x1, x2 = get_line(coord, i+I, j+J, nx+1, ny+1)
        return (x1 = x1, x2 = x2, equal_points = false)
    end

    function interpolate_line(I, J, L)
        pl = pillar_line(I, J)
        return interp_coord(pl, L)
    end

    l_11, t_11 = get_pair(0, 0)
    l_12, t_12 = get_pair(0, 1)
    l_21, t_21 = get_pair(1, 0)
    l_22, t_22 = get_pair(1, 1)

    pt_11 = interpolate_line(0, 0, l_11)
    pt_12 = interpolate_line(0, 1, l_12)
    pt_21 = interpolate_line(1, 0, l_21)
    pt_22 = interpolate_line(1, 1, l_22)


    A_1 = norm(cross(pt_21 - pt_11, pt_12 - pt_11), 2)
    A_2 = norm(cross(pt_21 - pt_22, pt_12 - pt_22), 2)
    area = (A_1 + A_2)/2.0

    d_11 = t_11 - l_11
    d_12 = t_12 - l_12
    d_21 = t_21 - l_21
    d_22 = t_22 - l_22

    d_avg = 0.25*(d_11 + d_12 + d_21 + d_22)
    return d_avg*area
end

function parse_physics_types(datafile)
    runspec = datafile["RUNSPEC"]
    props = datafile["PROPS"]
    has(name) = haskey(runspec, name) && runspec[name]
    has_wat = has("WATER")
    has_oil = has("OIL")
    has_gas = has("GAS")
    has_disgas = has("DISGAS")
    has_vapoil = has("VAPOIL")

    is_immiscible = !has_disgas && !has_vapoil
    is_compositional = haskey(runspec, "COMPS")

    phases = []
    rhoS = Vector{Float64}()
    pvt = []
    if haskey(props, "DENSITY")
        deck_density = props["DENSITY"]
        if deck_density isa Matrix
            deck_density = JutulDarcy.flat_region_expand(deck_density, 3)
        end
        if length(deck_density) > 1
            @warn "Multiple PVT regions found. Picking first one." deck_density
        end
        deck_density = deck_density[1]
        rhoOS = deck_density[1]
        rhoWS = deck_density[2]
        rhoGS = deck_density[3]
    else
        @assert is_compositional
        rhoOS = rhoWS = rhoGS = 1.0
        has_oil = true
        has_gas = true
    end

    if has_wat
        push!(pvt, JutulDarcy.deck_pvt_water(props))
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

        cnames = props["CNAMES"]
        acf = props["ACF"]
        mw = props["MW"]
        p_c = props["PCRIT"]
        V_c = props["VCRIT"]
        T_c = props["TCRIT"]

        if haskey(props, "BIC")
            A_ij = props["BIC"]
        else
            A_ij = nothing
        end
        mp = MolecularProperty.(mw, p_c, T_c, V_c, acf)
        mixture = MultiComponentMixture(mp, A_ij = A_ij, names = cnames)
        if haskey(props, "EOS")
            eos_str = lowercase(props["EOS"])
            if eos_str == "pr"
                eos_type = PengRobinson()
            elseif eos_str == "srk"
                eos_type = SoaveRedlichKwong()
            elseif eos_str == "rk"
                eos_type = RedlichKwong()
            else
                @assert eos_str == "zj"
                eos_type = ZudkevitchJoffe()
            end
        else
            eos_type = PengRobinson()
        end
        eos = GenericCubicEOS(mixture, eos_type)
        sys = MultiPhaseCompositionalSystemLV(eos, phases, reference_densities = rhoS)
    else
        if has_oil
            push!(pvt, JutulDarcy.deck_pvt_oil(props))
            push!(phases, LiquidPhase())
            push!(rhoS, rhoOS)
        end

        if has_gas
            push!(pvt, JutulDarcy.deck_pvt_gas(props))
            push!(phases, VaporPhase())
            push!(rhoS, rhoGS)
        end
        if is_immiscible
            sys = ImmiscibleSystem(phases, reference_densities = rhoS)
        else
            has_water = length(pvt) == 3
            oil_pvt = pvt[1 + has_water]
            if oil_pvt isa JutulDarcy.PVTO
                rs_max = JutulDarcy.saturated_table(oil_pvt)
            else
                rs_max = nothing
            end
            gas_pvt = pvt[2 + has_water]
            if gas_pvt isa JutulDarcy.PVTG
                rv_max = JutulDarcy.saturated_table(gas_pvt)
            else
                rv_max = nothing
            end
            sys = JutulDarcy.StandardBlackOilSystem(rs_max = rs_max, rv_max = rv_max, phases = phases,
                                        reference_densities = rhoS)
        end
    end

    return (system = sys, pvt = pvt)
end

function parse_schedule(domain, sys, datafile; kwarg...)
    schedule = datafile["SCHEDULE"]
    props = datafile["PROPS"]
    return parse_schedule(domain, datafile["RUNSPEC"], props, schedule, sys; kwarg...)
end

function parse_control_steps(runspec, props, schedule, sys)
    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)

    wells = schedule["WELSPECS"]
    steps = schedule["STEPS"]
    if length(keys(steps[end])) == 0
        # Prune empty final record
        steps = steps[1:end-1]
    end


    tstep = Vector{Float64}()
    cstep = Vector{Int}()
    compdat = Dict{String, OrderedDict}()
    controls = Dict{String, Any}()
    limits = Dict{String, Any}()
    streams = Dict{String, Any}()
    well_injection = Dict{String, Any}()
    for k in keys(wells)
        compdat[k] = OrderedDict{NTuple{3, Int}, Any}()
        controls[k] = DisabledControl()
        limits[k] = nothing
        streams[k] = nothing
        well_injection[k] = nothing
    end
    all_compdat = []
    all_controls = []
    all_limits = []

    if haskey(runspec, "START")
        start_date = runspec["START"]
    else
        start_date = missing
    end
    current_time = 0.0
    function add_dt!(dt, ctrl_ix)
        @assert dt > 0.0
        push!(tstep, dt)
        push!(cstep, ctrl_ix)
    end

    skip = ("WEFAC", "WELTARG", "WELLSTRE", "WINJGAS", "GINJGAS", "GRUPINJE", "WELLINJE")
    for (ctrl_ix, step) in enumerate(steps)
        found_time = false
        streams = parse_well_streams_for_step(step, props)
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
                        entry[(I, J, K)] = (
                            open = flag == "OPEN",
                            satnum = satnum,
                            WI = WI,
                            diameter = diam,
                            Kh = Kh,
                            skin = skin,
                            dir = dir,
                            ctrl = ctrl_ix)
                    end
                end
            elseif key in ("WCONINJE", "WCONPROD", "WCONHIST", "WCONINJ", "WCONINJH")
                for wk in kword
                    name = wk[1]
                    controls[name], limits[name] = keyword_to_control(sys, streams, wk, key)
                end
            elseif key in ("WEFAC", "WELTARG")
                @warn "$key Not supported properly."
            elseif key in skip
                # Already handled
            else
                @warn "Unhandled well keyword $key"
            end
        end
        push!(all_compdat, deepcopy(compdat))
        push!(all_controls, deepcopy(controls))
        push!(all_limits, deepcopy(limits))
        if !found_time
            error("Did not find supported time kw in step $ctrl_ix: Keys were $(keys(step))")
        end
    end
    return (dt = tstep, control_step = cstep, controls = all_controls, completions = all_compdat, limits = all_limits)
end

function parse_well_streams_for_step(step, props)
    if haskey(props, "STCOND")
        std = props["STCOND"]
        T = convert_to_si(std[1], :Celsius)
        p = std[2]
    else
        @debug "Defaulted STCOND..."
        T = convert_to_si(15.56, :Celsius)
        p = convert_to_si(1.0, :atm)
    end

    streams = Dict{String, Any}()
    well_streams = Dict{String, String}()
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

function keyword_to_control(sys, streams, kw, k::String)
    return keyword_to_control(sys, streams, kw, Val(Symbol(k)))
end

function keyword_to_control(sys, streams, kw, ::Val{:WCONPROD})
    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)

    flag = kw[2]
    ctrl = kw[3]
    orat = kw[4]
    wrat = kw[5]
    grat = kw[6]
    lrat = kw[7]
    bhp = kw[9]
    return producer_control(sys, flag, ctrl, orat, wrat, grat, lrat, bhp)
end

function keyword_to_control(sys, streams, kw, ::Val{:WCONHIST})
    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)
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

    return producer_control(sys, flag, ctrl, orat, wrat, grat, lrat, bhp, is_hist = true)
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

function producer_control(sys, flag, ctrl, orat, wrat, grat, lrat, bhp; is_hist = false)
    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)

    if flag == "SHUT" || flag == "STOP"
        ctrl = DisabledControl()
        lims = nothing
    else
        @assert flag == "OPEN"
        if ctrl == "LRAT"
            t = SurfaceLiquidRateTarget(-lrat)
        elseif ctrl == "WRAT"
            t = SurfaceWaterRateTarget(-wrat)
        elseif ctrl == "ORAT"
            t = SurfaceOilRateTarget(-orat)
        elseif ctrl == "GRAT"
            t = SurfaceGasRateTarget(-grat)
        elseif ctrl == "BHP"
            t = BottomHolePressureTarget(bhp)
        elseif ctrl == "RESV"
            t = -(wrat + orat + grat)
            if is_hist
                ctrl = HistoricalReservoirVoidageTarget(t)
            else
                ctrl = ReservoirVoidageTarget(t)
            end
        else
            error("$ctype control not supported")
        end
        ctrl = ProducerControl(t)
        lims = producer_limits(bhp = bhp, orat = orat, wrat = wrat, grat = grat, lrat = lrat)
    end
    return (ctrl, lims)
end

function injector_limits(; bhp = Inf, surface_rate = Inf, reservoir_rate = Inf)
    lims = Dict{Symbol, Any}()
    if isfinite(bhp)
        lims[:bhp] = bhp
    end
    if isfinite(surface_rate)
        lims[:rate] = surface_rate
    end
    if !isinf(reservoir_rate)
        @warn "Non-defaulted reservoir rate limit not supported: $reservoir_rate"
    end
    return NamedTuple(pairs(lims))
end

function injector_control(sys, streams, name, flag, type, ctype, surf_rate, res_rate, bhp; is_hist = false)
    if flag == "SHUT" || flag == "STOP"
        ctrl = DisabledControl()
        lims = nothing
    else
        @assert flag == "OPEN"
        if ctype == "RATE"
            t = TotalRateTarget(surf_rate)
        elseif ctype == "BHP"
            t = BottomHolePressureTarget(bhp)
        else
            # RESV, GRUP, THP
            error("$ctype control not supported")
        end
        rho, mix = select_injector_mixture_spec(sys, name, streams, type)
        ctrl = InjectorControl(t, mix, density = rho)
        lims = injector_limits(bhp = bhp, surface_rate = surf_rate, reservoir_rate = res_rate)
    end
    return (ctrl, lims)
end

function select_injector_mixture_spec(sys::Union{ImmiscibleSystem, StandardBlackOilSystem}, name, streams, type)
    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)
    mix = Float64[]
    rho = 0.0
    for (phase, rho_ph) in zip(phases, rho_s)
        if phase == LiquidPhase()
            v = Float64(type == "OIL")
        elseif phase == AqueousPhase()
            v = Float64(type == "WATER")
        else
            @assert phase isa VaporPhase
            v = Float64(type == "GAS")
        end
        rho += rho_ph*v
        push!(mix, v)
    end
    @assert sum(mix) ≈ 1.0
    return (rho, mix)
end

function select_injector_mixture_spec(sys::CompositionalSystem, name, streams, type)
    eos = sys.equation_of_state
    props = eos.mixture.properties
    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)
    mix = Float64[]
    rho = 0.0
    # Well stream will be molar fracitons.
    offset = Int(has_other_phase(sys))
    ncomp = number_of_components(sys)
    mix = zeros(Float64, ncomp)
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
    for i in 1:ncomp
        mix[i+offset] = z_mass[i]
    end
    @assert sum(mix) ≈ 1.0 "Sum of mixture was $(sum(mix)) != 1 for mole mixture $(z) as mass $z_mass"

    flash_cond = (p = cond.p, T = cond.T, z = z)
    flash = MultiComponentFlash.flashed_mixture_2ph(eos, flash_cond)
    rho_l, rho_v = MultiComponentFlash.mass_densities(eos, cond.p, cond.T, flash)
    S_l, S_v = MultiComponentFlash.phase_saturations(eos, cond.p, cond.T, flash)
    rho = S_l*rho_l + S_v*rho_v
    return (rho, mix)
end

function keyword_to_control(sys, streams, kw, ::Val{:WCONINJE})
    # TODO: Expand to handle mixture etc.
    name = kw[1]
    type = kw[2]
    flag = kw[3]
    ctype = kw[4]
    surf_rate = kw[5]
    res_rate = kw[6]
    bhp = kw[7]
    return injector_control(sys, streams, name, flag, type, ctype, surf_rate, res_rate, bhp)
end

function keyword_to_control(sys, streams, kw, ::Val{:WCONINJH})
    name = kw[1]
    type = kw[2]
    flag = kw[3]
    surf_rate = kw[4]
    bhp = kw[5]
    ctype = kw[12]
    # TODO: Expand to handle mixture etc.
    res_rate = Inf
    return injector_control(sys, streams, name, flag, type, ctype, surf_rate, res_rate, bhp, is_hist = true)
end
