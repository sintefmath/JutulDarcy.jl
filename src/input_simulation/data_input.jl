export simulate_data_file, case_from_data_input

"""
    simulate_data_file(inp; parse_arg = NamedTuple(), kwarg...)

Simulate standard input file (with extension .DATA). `inp` can either be the
output from `case_from_data_input` or a String for the path of an input file.

Additional arguments are passed onto `simulate_reservoir`. Extra inputs to the
parser can be sent as a `parse_arg` `NamedTuple`.
"""
function simulate_data_file(fn::String; parse_arg = NamedTuple(), kwarg...)
    data = parse_data_file(fn; parse_arg...)
    return simulate_data_file(data; kwarg...)
end

function simulate_data_file(data::AbstractDict; setup_arg = NamedTuple(), kwarg...)
    case = case_from_data_input(data; setup_arg...)
    return simulate_reservoir(case; kwarg...)
end

function case_from_data_input(datafile; kwarg...)
    sys, pvt = parse_physics_types(datafile)
    domain = parse_reservoir(datafile)
    wells, controls, limits, cstep, dt, well_mul = parse_schedule(domain, sys, datafile)

    model, parameters0 = setup_reservoir_model(domain, sys; wells = wells, kwarg...)

    for (k, submodel) in pairs(model.models)
        if submodel.system isa MultiPhaseSystem
            # Modify secondary variables
            svar = submodel.secondary_variables
            # PVT
            pvt = tuple(pvt...)
            rho = DeckDensity(pvt)
            if sys isa StandardBlackOilSystem
                set_secondary_variables!(submodel, ShrinkageFactors = JutulDarcy.DeckShrinkageFactors(pvt))
            end
            mu = DeckViscosity(pvt)
            set_secondary_variables!(submodel, PhaseViscosities = mu, PhaseMassDensities = rho)
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
    forces = parse_forces(model, wells, controls, limits, cstep, dt, well_mul)
    state0 = parse_state0(model, datafile)
    return JutulCase(model, dt, forces, state0 = state0, parameters = parameters)
end

function parse_well_from_compdat(domain, wname, v, wspecs)
    wc, WI, open = compdat_to_connection_factors(domain, v)
    rd = wspecs.ref_depth
    if isnan(rd)
        rd = nothing
    end
    W = setup_well(domain, wc, name = Symbol(wname), WI = WI, reference_depth = rd)
    return (W, WI, open)
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

function parse_schedule(domain, runspec, schedule, sys)
    dt, cstep, controls, completions, limits = parse_control_steps(runspec, schedule, sys)
    ncomp = length(completions)
    wells = []
    well_mul = []
    for (k, v) in pairs(completions[end])
        W, WI, open = parse_well_from_compdat(domain, k, v, schedule["WELSPECS"][k])
        wi_mul = zeros(length(WI), ncomp)
        for (i, c) in enumerate(completions)
            compdat = c[k]
            _, WI, open = compdat_to_connection_factors(domain, compdat)
            @. wi_mul[:, i] = (WI*open)/WI
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

    nph = number_of_phases(sys)
    sat = zeros(nph, nc)
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

    if is_blackoil
        if ismissing(sw)
            sw = zeros(nc)
        end
        @assert !ismissing(sg)
        if nph == 3
            init[:ImmiscibleSaturation] = sw
        end
        F_rs = sys.rs_max
        F_rv = sys.rv_max

        so = s_rem
        init[:BlackOilUnknown] = map(
            (w,  o,   g, r,  v, p) -> JutulDarcy.blackoil_unknown_init(F_rs, F_rv, w, o, g, r, v, p),
             sw, so, sg, rs, rv, pressure)
    else
        error("Not implemented yet.")
    end
    return init
end

function parse_reservoir(data_file)
    grid = data_file["GRID"]
    coord = grid["COORD"]
    zcorn = grid["ZCORN"]
    actnum = grid["ACTNUM"]
    cartdims = grid["cartDims"]
    primitives = JutulDarcy.cpgrid_primitives(coord, zcorn, cartdims, actnum = actnum)
    G = JutulDarcy.grid_from_primitives(primitives)

    nc = number_of_cells(G)
    ix = G.cell_map
    # TODO: PERMYY etc for full tensor perm
    perm = zeros(3, nc)
    poro = zeros(nc)
    for (i, c) in enumerate(G.cell_map)
        perm[1, i] = grid["PERMX"][c]
        perm[2, i] = grid["PERMY"][c]
        perm[3, i] = grid["PERMZ"][c]
        poro[i] = grid["PORO"][c]
    end
    satnum = JutulDarcy.InputParser.table_region(data_file, :saturation, active = G.cell_map)
    eqlnum = JutulDarcy.InputParser.table_region(data_file, :equil, active = G.cell_map)

    domain = reservoir_domain(G, permeability = perm, porosity = poro, satnum = satnum, eqlnum = eqlnum)
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
    # TODO: compositional detection
    # is_compositional = haskey(mrst_data, "mixture")

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

    return (system = sys, pvt = pvt)
end

function parse_schedule(domain, sys, datafile)
    schedule = datafile["SCHEDULE"]
    return parse_schedule(domain, datafile["RUNSPEC"], schedule, sys)
end

function parse_control_steps(runspec, schedule, sys)
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
    for k in keys(wells)
        compdat[k] = OrderedDict{NTuple{3, Int}, Any}()
        controls[k] = nothing
        limits[k] = nothing
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

    for (ctrl_ix, step) in enumerate(steps)
        found_time = false
        for (key, kword) in pairs(step)
            if key == "DATES"
                @assert !ismissing(start_date)
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
            elseif key in ("WCONINJE", "WCONPROD", "WCONHIST", "WCONINJ")
                for wk in kword
                    name = wk[1]
                    controls[name], limits[name] = keyword_to_control(sys, wk, key)
                end
            elseif key in ("WEFAC", "WELTARG")
                @warn "$key Not supported properly."
            else
                error("Unhandled keyword $key")
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


function keyword_to_control(sys, kw, k::String)
    return keyword_to_control(sys, kw, Val(Symbol(k)))
end

function keyword_to_control(sys, kw, ::Val{:WCONPROD})
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

function keyword_to_control(sys, kw, ::Val{:WCONHIST})
    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)
    # 1 name
    error("Not implemented yet.")
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

function producer_control(sys::Union{ImmiscibleSystem, StandardBlackOilSystem}, flag, ctrl, orat, wrat, grat, lrat, bhp)
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

function injector_control(sys::Union{ImmiscibleSystem, StandardBlackOilSystem}, flag, type, ctype, surf_rate, res_rate, bhp)
    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)

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
        ctrl = InjectorControl(t, mix, density = rho)
        lims = injector_limits(bhp = bhp, surface_rate = surf_rate, reservoir_rate = res_rate)
    end
    return (ctrl, lims)
end

function keyword_to_control(sys, kw, ::Val{:WCONINJE})
    type = kw[2]
    flag = kw[3]
    ctype = kw[4]
    surf_rate = kw[5]
    res_rate = kw[6]
    bhp = kw[7]
    return injector_control(sys, flag, type, ctype, surf_rate, res_rate, bhp)
end
