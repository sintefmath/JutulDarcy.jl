
"""
    equilibriate_state(model, contacts)

Equilibrates the state of the given model based on the provided contacts.

# Arguments
- `model`: The model whose state needs to be equilibrated.
- `contacts`: The nph contact depths.

# Keyword Arguments
- `datum_depth`: The reference depth for the datum.
- `datum_pressure`: The pressure at the datum depth.
- `cells`: The cells to be equilibrated.
- `rs`: Solution gas-oil ratio (blackoil).
- `rv`: Vapor-oil ratio (blackoil).
- `composition`: The composition vs depth (compositional).
- `kwarg`: Additional keyword arguments.

# Returns
- The equilibrated state of the model.
"""
function equilibriate_state(model, contacts,
        datum_depth = missing,
        datum_pressure = JutulDarcy.DEFAULT_MINIMUM_PRESSURE;
        cells = missing,
        rs = missing,
        rv = missing,
        composition = missing,
        kwarg...
    )
    model = reservoir_model(model)
    D = model.data_domain
    G = physical_representation(D)
    cc = D[:cell_centroids][3, :]
    if ismissing(cells)
        cells = 1:number_of_cells(G)
        pts = cc
    else
        pts = view(cc, cells)
    end

    if ismissing(datum_depth)
        datum_depth = minimum(pts)
    end
    sys = flow_system(model.system)

    init = Dict{Symbol, Any}()
    init = equilibriate_state!(init, pts, model, sys, contacts, datum_depth, datum_pressure;
        cells = cells, rv = rv, rs = rs, composition = composition, kwarg...
    )

    is_blackoil = sys isa StandardBlackOilSystem
    if is_blackoil
        sat = init[:Saturations]
        pressure = init[:Pressure]
        nph, nc = size(sat)

        if has_disgas(sys)
            if ismissing(rs)
                Rs = zeros(nc)
            else
                Rs = rs.(pts)
            end
            init[:Rs] = Rs
        end
        if has_vapoil(sys)
            if ismissing(rv)
                Rv = zeros(nc)
            else
                Rv = rv.(pts)
            end
            init[:Rv] = Rv
        end
    end
    is_compositional = sys isa CompositionalSystem
    if is_compositional
        @assert !ismissing(composition)
        has_wat = has_other_phase(sys)
        nc = length(init[:Pressure])
        ncomp = number_of_components(sys) - has_wat
        zmf = zeros(ncomp, nc)
        for i in 1:nc
            zmf[:, i] = composition(pts[i])
        end
        init[:OverallMoleFractions] = zmf
    end
    if (is_compositional || is_blackoil) && has_other_phase(sys)
        init[:ImmiscibleSaturation] = init[:Saturations][1, :]
    end
    return init
end

function equilibriate_state(model, equil::EquilibriumRegion)
    sys = model.system
    phases = get_phases(sys)
    phase_ix = [i for i in phase_indices(sys)]
    nph = length(phases)
    nph in (1, 2, 3) || error("Only 1, 2, or 3 phases are supported for equilibriation with EquilibriumRegion.")

    has_water = AqueousPhase() in phases
    has_oil = LiquidPhase() in phases
    has_gas = VaporPhase() in phases

    function check_pair(c, pc_c, name)
        !isnan(c) || throw(ArgumentError("$name cannot be defaulted for this system."))
        isfinite(pc_c) || throw(ArgumentError("Capillary pressure at contact (pc_$name) is not finite."))
    end

    contacts = Float64[]
    contacts_pc = Float64[]
    if has_water && has_oil
        push!(contacts, equil.woc)
        push!(contacts_pc, equil.pc_woc)
        check_pair(equil.woc, equil.pc_woc, "woc")
    elseif has_water && has_gas
        push!(contacts, equil.wgc)
        push!(contacts_pc, equil.pc_wgc)
        check_pair(equil.wgc, equil.pc_wgc, "wgc")
    end
    if has_oil && has_gas
        goc = equil.goc
        push!(contacts, equil.goc)
        push!(contacts_pc, equil.pc_goc)
        check_pair(equil.goc, equil.pc_goc, "goc")
    end
    model = reservoir_model(model)
    init = equilibriate_state(model,
        contacts,
        equil.datum_depth,
        equil.datum_pressure;
        contacts_pc = contacts_pc,
        T_z = equil.temperature_vs_depth,
        rs = equil.rs_vs_depth,
        rv = equil.rv_vs_depth,
        composition = equil.composition_vs_depth,
        cells = equil.cells,
        density_function = equil.density_function,
        pvtnum = equil.pvtnum,
        satnum = equil.satnum,
        equil.kwarg...
    )
    init[:Saturations] = init[:Saturations][phase_ix, :]
    return init
end

function equilibriate_state!(init, depths, model, sys, contacts, depth, datum_pressure;
        cells = 1:length(depths),
        rs = missing,
        rv = missing,
        pc = missing,
        composition = missing,
        T_z = missing,
        s_min = missing,
        contacts_pc = missing,
        pvtnum = 1,
        satnum = 1,
        sw = missing,
        output_pressures = false,
        density_function = missing,
        kwarg...
    )
    if ismissing(pc)
        pc_def = get(model.secondary_variables, :CapillaryPressure, nothing)
        if !isnothing(pc_def)
            pc = []
            for (i, f) in enumerate(pc_def.pc)
                f = table_by_region(f, satnum)
                s = copy(f.X)
                cap = copy(f.F)
                if s[1] < 0
                    s = s[2:end]
                    cap = cap[2:end]
                end
                ix = unique(i -> cap[i], 1:length(cap))
                s = s[ix]
                cap = cap[ix]
                if length(s) == 1
                    push!(s, s[end])
                    push!(cap, cap[end]+1.0)
                end
                push!(pc, (s = s, pc = cap))
            end
        end
    end
    if ismissing(contacts_pc)
        contacts_pc = zeros(number_of_phases(sys)-1)
    end
    zmin = minimum(depths) - 1.0
    zmax = maximum(depths) + 1.0

    nph = number_of_phases(sys)

    @assert length(contacts) == nph-1
    rho = model.secondary_variables[:PhaseMassDensities]
    if rho isa Pair
        rho = last(rho)
    end

    reg = Int[pvtnum]
    # Set up a mock state for evaluation
    if haskey(model.data_domain, :pvtnum)
        fake_cell_ix = [findfirst(isequal(pvtnum), model.data_domain[:pvtnum])]
    else
        fake_cell_ix = [1]
    end
    fake_state = JutulStorage()
    fake_state[:Pressure] = [NaN]
    fake_state[:PhaseMassDensities] = zeros(nph, 1)
    fake_state[:Temperature] = [NaN]

    density_f(p, z, ph) = equilibrium_phase_density(p, z, ph, rho, T_z, fake_state, model, rs, rv, composition, fake_cell_ix, reg)
    # Find the reference phase. It is either liquid
    ref_phase = get_reference_phase_index(model.system)
    if ismissing(density_function)
        density_function = density_f
    end
    pressures = determine_hydrostatic_pressures(depths, depth, zmin, zmax, contacts, datum_pressure, density_function, contacts_pc, ref_phase)
    if nph > 1
        relperm = model.secondary_variables[:RelativePermeabilities]
        if relperm isa Pair
            relperm = last(relperm)
        end
        s, pc = determine_saturations(depths, contacts, pressures; pc = pc, s_min = s_min, kwarg...)
        if !ismissing(sw)
            nph = size(s, 1)
            for i in axes(s, 2)
                if pressures[2, i] - pressures[1, i] < 0.0
                    continue
                end
                s0 = s[:, i]
                sw_i = sw[i]
                s[1, i] = sw_i
                s_rem = 0.0
                for ph in 2:nph
                    s_rem += s[ph, i]
                end
                if s_rem ≈ 0
                    s[2:end, i] .= 0.0
                else
                    scale = (1 - sw_i)/max(s_rem, 1e-20)
                    for ph in 2:nph
                        s[ph, i] *= scale
                    end
                    sT = sum(s[:, i])
                    if !isapprox(sT, 1.0, atol = 1e-8)
                        @warn "$s0 $(sw[i]) $scale Was $sT, not 1.0"
                    end
                end
            end
        end
        if relperm isa ReservoirRelativePermeabilities
            nc_total = number_of_cells(model.domain)
            T_S = eltype(s)
            kr = zeros(T_S, nph, nc_total)
            s_eval = zeros(T_S, nph, nc_total)
            s_eval[:, cells] .= s
            phases = get_phases(sys)
            phase_ind = phase_indices(model.system)
            if length(phases) == 3
                swcon = zeros(nc_total)
                if AqueousPhase() in phases && !ismissing(s_min)
                    swcon[cells] .= s_min[1]
                end
                for c in cells
                    update_three_phase_relperm!(kr, relperm, phase_ind, s_eval, nothing, c, swcon[c], nothing, nothing)
                end
            else
                ph = relperm.phases
                if ph == :wo
                    kr1, kr2 = relperm.krw, relperm.krow
                elseif ph == :og
                    kr1, kr2 = relperm.krog, relperm.krg
                else
                    @assert ph == :wg
                    kr1, kr2 = relperm.krw, relperm.krg
                end
                for c in cells
                    update_two_phase_relperm!(kr, relperm, kr1, kr2, nothing, nothing, phase_ind, s_eval, nothing, c, nothing, nothing)
                end
            end
            kr = kr[:, cells]
        else
            jutul_message("Initialization", "Rel. perm. is not a ReservoirRelativePermeabilities, skipping mobility check for reference pressure.")
            kr = copy(s)
        end
        init[:Saturations] = s
        init[:Pressure] = init_reference_pressure(pressures, contacts, kr, pc, ref_phase)
    else
        p = copy(vec(pressures))
        init[:Pressure] = p
        init[:Saturations] = ones(1, length(p))
    end
    if output_pressures
        init[:EquilibriationPressures] = pressures
    end

    if !ismissing(T_z)
        init[:Temperature] = T_z.(depths)
    end

    return init
end

function equilibrium_phase_density(p, z, ph, rho, T_z, fake_state, model, rs, rv, composition, fake_cell_ix, reg)
    pvtnum = only(reg)
    sys = model.system
    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)
    disgas = JutulDarcy.has_disgas(sys)
    vapoil = JutulDarcy.has_vapoil(sys)

    if disgas || vapoil
        if JutulDarcy.has_other_phase(sys)
            _, rhoOS, rhoGS = rho_s
        else
            rhoOS, rhoGS = rho_s
        end
    end
    if rho isa ThreePhaseCompositionalDensitiesLV || rho isa TwoPhaseCompositionalDensities
        if phases[ph] == AqueousPhase()
            phase_density = rho_s[ph]*JutulDarcy.shrinkage(rho.immiscible_pvt, p)
        else
            @assert !ismissing(composition) "Composition must be present for equilibrium calculations for compositional models."
            @assert !ismissing(T_z) "Temperature must be present for equilibrium calculations for compositional models."
            z_i = Vector{Float64}(composition(z))
            T = T_z(z)
            eos = model.system.equation_of_state
            f = MultiComponentFlash.flashed_mixture_2ph(eos, (p = p, T = T, z = z_i))
            rho_l, rho_v = mass_densities(eos, p, T, f)
            if phases[ph] == VaporPhase() || rho_l == 0
                phase_density = rho_v
            else
                phase_density = rho_l
            end
        end
    elseif rho isa BrineCO2MixingDensities
        T = T_z(z)
        phase_density = rho.tab(p, T)[ph]
    elseif rho isa DeckPhaseMassDensities
        pvt = rho.pvt
        pvt_i = pvt[ph]
        if phases[ph] == LiquidPhase() && disgas
            rs_max = table_by_region(sys.rs_max, pvtnum)
            Rs = min(rs(z), rs_max(p))
            b = JutulDarcy.shrinkage(pvt_i, reg, p, Rs, 1)
            phase_density = b*(rhoOS + Rs*rhoGS)
        elseif phases[ph] == VaporPhase() && vapoil
            rv_max = table_by_region(sys.rv_max, pvtnum)
            Rv = min(rv(z), rv_max(p))
            b = JutulDarcy.shrinkage(pvt_i, reg, p, Rv, 1)
            phase_density = b*(rhoGS + Rv*rhoOS)
        else
            phase_density = rho_s[ph]*JutulDarcy.shrinkage(pvt_i, reg, p, 1)
        end
    else
        rho_val = fake_state[:PhaseMassDensities]
        fake_state[:Pressure][1] = p
        if !ismissing(T_z)
            fake_state[:Temperature][1] = T_z(z)
        end
        Jutul.update_secondary_variable!(rho_val, rho, model, fake_state, fake_cell_ix)
        phase_density = rho_val[ph]
    end
    return phase_density
end

function parse_state0_equil(model, datafile; normalize = :sum)
    sys = model.system
    d = model.data_domain

    has_water = haskey(datafile["RUNSPEC"], "WATER")
    if sys isa CompositionalSystem
        has_oil = has_gas = true
    else
        has_oil = haskey(datafile["RUNSPEC"], "OIL")
        has_gas = haskey(datafile["RUNSPEC"], "GAS")
    end

    is_co2 = haskey(datafile["RUNSPEC"], "JUTUL_CO2BRINE")
    is_single_phase = (has_water + has_oil + has_gas) == 1

    has_sat_reg = haskey(d, :satnum)
    ncells = number_of_cells(d)
    if has_sat_reg
        satnum = d[:satnum]
    else
        satnum = ones(Int, ncells)
    end
    if haskey(d, :pvtnum)
        pvtnum = d[:pvtnum]
    else
        pvtnum = ones(Int, ncells)
    end
    if haskey(d, :eqlnum)
        eqlnum = d[:eqlnum]
    else
        eqlnum = ones(Int, ncells)
    end

    if is_single_phase
        has_pc = false
        pc_functions = missing
        krw_fn = missing
    else
        has_pc = haskey(model.secondary_variables, :CapillaryPressure)
        if has_pc
            pcvar = model.secondary_variables[:CapillaryPressure]
            pc_functions = pcvar.pc
            if has_sat_reg
                @assert !isnothing(pcvar.regions)
                @assert pcvar.regions == satnum
            end
        else
            pc_functions = missing
        end

        if has_water
            kr_var = model.secondary_variables[:RelativePermeabilities]
            if kr_var isa ReservoirRelativePermeabilities
                krw_fn = kr_var.krw
                if has_sat_reg
                    @assert !isnothing(kr_var.regions)
                    @assert kr_var.regions == satnum
                end
            else
                krw_fn = missing
            end
        else
            krw_fn = missing
        end
    end

    init = Dict{Symbol, Any}()
    sol = datafile["SOLUTION"]
    props = datafile["PROPS"]

    G = physical_representation(model.data_domain)
    nc = number_of_cells(G)
    nph = number_of_phases(model.system)
    actnum_ix = G.cell_map
    if isnothing(actnum_ix)
        actnum_ix = 1:nc
    end
    is_blackoil = sys isa StandardBlackOilSystem
    disgas = JutulDarcy.has_disgas(model.system)
    vapoil = JutulDarcy.has_vapoil(model.system)

    equil = sol["EQUIL"]
    nequil = GeoEnergyIO.InputParser.number_of_tables(datafile, :eqlnum)
    npvt = GeoEnergyIO.InputParser.number_of_tables(datafile, :pvtnum)
    nsat = GeoEnergyIO.InputParser.number_of_tables(datafile, :satnum)

    @assert length(equil) == nequil
    inits = []
    inits_cells = []
    for ereg in 1:nequil
        T_z = missing
        if haskey(sol, "RTEMP") || haskey(props, "RTEMP")
            if haskey(props, "RTEMP")
                rtmp = props["RTEMP"][1]
            else
                rtmp = sol["RTEMP"][1]
            end
            Ti = convert_to_si(rtmp, :Celsius)
            T_z = z -> Ti
        else
            for sect in [props, sol]
                if haskey(sect, "TEMPVD") || haskey(sect, "RTEMPVD")
                    if haskey(sect, "TEMPVD")
                        tvd_kw = sect["TEMPVD"][ereg]
                    else
                        tvd_kw = sect["RTEMPVD"][ereg]
                    end
                    z = vec(tvd_kw[:, 1])
                    Tvd = vec(tvd_kw[:, 2] .+ 273.15)
                    T_z = get_1d_interpolator(z, Tvd)
                end
            end
        end
        eq = equil[ereg]
        cells_eqlnum = findall(isequal(ereg), eqlnum)
        for sreg in 1:nsat
            cells_satnum = findall(isequal(sreg), satnum)
            cells_sat_and_pvt = intersect_sorted(cells_satnum, cells_eqlnum)
            for preg in 1:npvt
                cells_pvtnum = findall(isequal(preg), pvtnum)
                cells = intersect_sorted(cells_pvtnum, cells_sat_and_pvt)
                ncells_reg = length(cells)
                if ncells_reg == 0
                    continue
                end
                actnum_ix_for_reg = actnum_ix[cells]
                datum_depth = eq[1]
                datum_pressure = eq[2]

                woc = eq[3]
                woc_pc = eq[4]
                goc = eq[5]
                goc_pc = eq[6]
                rs_method = eq[7]
                rv_method = eq[8]
                # Contact depths
                s_max = 1.0
                s_min = 0.0

                non_connate = ones(ncells_reg)
                s_max = Vector{Float64}[]
                s_min = Vector{Float64}[]
                if has_pc
                    pc = []
                    for (i, f) in enumerate(pc_functions)
                        f = table_by_region(f, sreg)
                        s = copy(f.X)
                        cap = copy(f.F)
                        if s[1] < 0
                            s = s[2:end]
                            cap = cap[2:end]
                        end
                        ix = unique(i -> cap[i], 1:length(cap))
                        s = s[ix]
                        cap = cap[ix]
                        if length(s) == 1
                            push!(s, s[end])
                            push!(cap, cap[end]+1.0)
                        end
                        push!(pc, (s = s, pc = cap))
                    end
                else
                    pc = nothing
                end
                if is_single_phase
                    push!(s_min, zeros(ncells_reg))
                    push!(s_max, ones(ncells_reg))
                else
                    if has_water
                        if ismissing(krw_fn)
                            push!(s_min, zeros(ncells_reg))
                            push!(s_max, ones(ncells_reg))
                        else
                            krw = table_by_region(krw_fn, sreg)
                            if haskey(datafile, "PROPS") && haskey(datafile["PROPS"], "SWL")
                                swl = vec(datafile["PROPS"]["SWL"])
                                swcon = swl[actnum_ix_for_reg]
                            else
                                swcon = fill(krw.connate, ncells_reg)
                            end
                            push!(s_min, swcon)
                            push!(s_max, fill(krw.input_s_max, ncells_reg))
                            @. non_connate -= swcon
                        end
                    end
                    if has_oil
                        push!(s_min, zeros(ncells_reg))
                        push!(s_max, non_connate)
                    end
                    if has_gas
                        push!(s_min, zeros(ncells_reg))
                        push!(s_max, non_connate)
                    end
                end

                if nph == 1
                    contacts = Float64[]
                    contacts_pc = Float64[]
                elseif nph == 2
                    if is_co2
                        contacts = (woc, )
                        # TODO: Check sign here. Usually these models are
                        # initialized without CO2.
                        contacts_pc = (woc_pc, )
                    elseif has_oil && has_gas
                        contacts = (goc, )
                        contacts_pc = (goc_pc, )
                    else
                        contacts = (woc, )
                        contacts_pc = (woc_pc, )
                    end
                else
                    contacts = (woc, goc)
                    contacts_pc = (woc_pc, goc_pc)
                end

                if disgas || vapoil
                    rhoGS = map(x -> x[3], datafile["PROPS"]["DENSITY"])
                    rhoOS = map(x -> x[1], datafile["PROPS"]["DENSITY"])

                    Rs_scale = (rhoGS[preg]/rhoGS[1])*(rhoOS[preg]/rhoOS[1])
                    Rv_scale = 1.0/Rs_scale
                    if disgas
                        if rs_method <= 0
                            rs_max_at_contact = sys.rs_max[preg](datum_pressure)
                            rs = z -> rs_max_at_contact
                        else
                            if haskey(sol, "RSVD")
                                rsvd = sol["RSVD"][ereg]
                                z = rsvd[:, 1]
                                Rs = rsvd[:, 2]
                            else
                                @assert haskey(sol, "PBVD")
                                pbvd = sol["PBVD"][ereg]
                                z = pbvd[:, 1]
                                pb = vec(pbvd[:, 2])
                                Rs = sys.rs_max[preg].(pb)
                            end
                            rs = get_1d_interpolator(z, Rs_scale.*Rs)
                        end
                    else
                        rs = missing
                    end
                    if vapoil
                        if rv_method <= 0
                            rv_max_at_contact = sys.rv_max[preg](datum_pressure)
                            rv = z -> rv_max_at_contact
                        else
                            if haskey(sol, "PDVD")
                                @warn "PDVD not supported for RV initialization, setting to zero."
                                rv = z -> 0.0
                            elseif haskey(sol, "RVVD")
                                rvvd = sol["RVVD"][ereg]
                                z = rvvd[:, 1]
                                Rv = rvvd[:, 2]
                                rv = get_1d_interpolator(z, Rv_scale.*Rv)
                            else
                                # TODO: Should maybe be the pressure at GOC not datum
                                # depth, but that isn't computed at this stage.
                                rv_at_goc = Rv_scale*sys.rv_max[preg](datum_pressure)
                                rv = z -> rv_at_goc
                            end
                        end
                    else
                        rv = missing
                    end
                else
                    rs = rv = missing
                end
                if haskey(props, "ZMFVD")
                    ztab = props["ZMFVD"][ereg]
                    N = size(ztab, 2) - 1
                    Sz = SVector{N, Float64}
                    ztab_z = Float64[]
                    ztab_static = Sz[]
                    num_renormalized = 0
                    for row in 1:size(ztab, 1)
                        mf_i = Sz(ztab[row, 2:end])
                        if maximum(mf_i) > 1.0 || minimum(mf_i) < 0 || !(sum(mf_i) ≈ 1.0)
                            num_renormalized += 1
                        end
                        if normalize == :sum || normalize == true
                            mf_i = mf_i./sum(mf_i)
                        elseif normalize == :clampsum
                            mf_i = clamp.(mf_i, 0.0, 1.0)
                            mf_i = mf_i./sum(mf_i)
                        else
                            @assert normalize == false
                        end
                        z_i = ztab[row, 1]
                        push!(ztab_z, z_i)
                        push!(ztab_static, mf_i)
                        if row == 1
                            push!(ztab_z, z_i-1.0)
                            push!(ztab_static, mf_i)
                        end
                        if row == size(ztab, 1)
                            push!(ztab_z, z_i+1.0)
                            push!(ztab_static, mf_i)
                        end
                    end
                    if num_renormalized > 0 && normalize != false
                        jutul_message("Initialization", "$num_renormalized entries ZMFVD were normalized in equilibriation region $ereg.")
                    end
                    composition = get_1d_interpolator(ztab_z, ztab_static)
                else
                    composition = missing
                end
                if haskey(props, "SWATINIT")
                    sw = props["SWATINIT"][actnum_ix[cells]]
                    extra_arg = (sw = sw, )
                else
                    extra_arg = NamedTuple()
                end

                subinit = equilibriate_state(
                        model, contacts, datum_depth, datum_pressure;
                        composition = composition,
                        cells = cells,
                        pvtnum = preg,
                        contacts_pc = contacts_pc,
                        s_min = s_min,
                        s_max = s_max,
                        T_z = T_z,
                        rs = rs,
                        rv = rv,
                        pc = pc,
                        output_pressures = true,
                        extra_arg...
                    )
                push!(inits, subinit)
                push!(inits_cells, cells)
            end
        end
    end
    if length(inits) == 1
        init = only(inits)
    else
        # Handle multiple regions by merging each init
        init = Dict{Symbol, Any}()
        nc = number_of_cells(model.domain)
        touched = [false for i in 1:nc]
        for (k, v) in first(inits)
            T_v = eltype(v)
            if v isa AbstractVector
                init[k] = zeros(T_v, nc)
            else
                @assert v isa AbstractMatrix
                init[k] = zeros(T_v, size(v, 1), nc)
            end
        end
        for (subinit, cells) in zip(inits, inits_cells)
            for c in cells
                if touched[c]
                    @warn "Equils overlap for cell $c?"
                end
                touched[c] = true
            end

            for (k, v) in subinit
                fill_subinit!(init[k], cells, v)
            end
        end
        @assert all(touched) "Some cells are not initialized by equil: $(findall(!, touched))"
    end
    if haskey(props, "SWATINIT") && nph > 1 && haskey(model.secondary_variables, :CapillaryPressure)
        sw = props["SWATINIT"][actnum_ix]
        pc = model.secondary_variables[:CapillaryPressure]

        sat = init[:Saturations]
        pcval = zeros(eltype(sat), nph-1, nc)
        update_pc!(pcval, pc, model, sat, 1:nc)
        pressure_eql = init[:EquilibriationPressures]
        pc_scale = ones(eltype(sat), nph-1, nc)
        for i in 1:nc
            sw_i = sw[i]
            pc_actual = pcval[1, i]
            pw = pressure_eql[1, i]
            po = pressure_eql[2, i]
            pc_eql = pw - po
            if abs(pc_actual) ≈ 0.0 || po < pw
                continue
            end
            # p_o + pc_ow = p_w
            # -> p_w - p_o = pc_ow
            pc_scale[1, i] = pc_eql/pc_actual
        end
        model.secondary_variables[:CapillaryPressure] = ScaledCapillaryPressure(pc.pc, pc_scale, regions = pc.regions)
    end
    delete!(init, :EquilibriationPressures)
    return init
end

function intersect_sorted(a::Vector{T}, b::Vector{T}) where T
    na = length(a)
    nb = length(b)
    c = Vector{T}()
    ptr_a = ptr_b = 1
    while ptr_a <= na && ptr_b <= nb
        a_val = a[ptr_a]
        b_val = b[ptr_b]
        if a_val == b_val
            push!(c, a_val)
            ptr_a += 1
            ptr_b += 1
        elseif a_val < b_val
            ptr_a += 1
        else
            ptr_b += 1
        end
    end
    return c
end

function fill_subinit!(x::Vector, cells, v::Vector)
    @assert length(v) == length(cells)
    for (i, c) in enumerate(cells)
        x[c] = v[i]
    end
end

function fill_subinit!(x::Matrix, cells, v::Matrix)
    @assert size(x, 1) == size(v, 1)
    @assert size(v, 2) == length(cells)
    for (i, c) in enumerate(cells)
        for j in axes(x, 1)
            x[j, c] = v[j, i]
        end
    end
end

function init_reference_pressure(pressures, contacts, kr, pc, ref_ix = 2)
    nph, nc = size(kr)
    T = promote_type(eltype(pressures), eltype(contacts), eltype(kr))
    p = zeros(T, nc)
    ϵ = 1e-12
    for i in eachindex(p)
        p[i] = pressures[ref_ix, i]
        kr_ref = kr[ref_ix, i]
        @assert kr_ref >= -ϵ "Evaluated rel. perm. was $kr_ref for phase reference phase (index $ref_ix)."
        if kr[ref_ix, i] <= ϵ
            for ph in 1:nph
                if kr[ph, i] > ϵ
                    p[i] = pressures[ph, i]
                end
            end
        end
    end
    return p
end

function determine_hydrostatic_pressures(depths, depth, zmin, zmax, contacts, datum_pressure, density_f, contacts_pc, ref_ix = 0)
    T = promote_type(eltype(depths), typeof(depth), typeof(zmin), typeof(zmax), eltype(contacts), typeof(datum_pressure))
    nc = length(depths)
    nph = length(contacts) + 1
    if ref_ix == 0
        ref_ix = min(2, nph)
    end
    I_ref = phase_pressure_depth_table(depth, zmin, zmax, datum_pressure, density_f, ref_ix)
    pressures = zeros(T, nph, nc)
    pos = 1
    for ph in 1:nph
        if ph == ref_ix
            I = I_ref
        else
            contact = contacts[pos]
            datum_pressure_ph = I_ref(contact) + contacts_pc[pos]
            I = phase_pressure_depth_table(contact, zmin, zmax, datum_pressure_ph, density_f, ph)
            pos += 1
        end
        for (i, z) in enumerate(depths)
            pressures[ph, i] = I(z)
        end
    end
    return pressures
end



function phase_pressure_depth_table(depth, zmin, zmax, datum_pressure, density_f, phase)
    if zmin > depth
        z_up = Float64[]
        p_up = Float64[]
    else
        z_up, p_up = integrate_phase_density(depth, zmin, datum_pressure, density_f, phase)
        # Remove top point
        z_up = reverse!(z_up[2:end])
        p_up = reverse!(p_up[2:end])
    end
    if zmax < depth
        z_down = Float64[]
        p_down = Float64[]
    else
        z_down, p_down = integrate_phase_density(depth, zmax, datum_pressure, density_f, phase)
    end
    z = vcat(z_up, z_down)
    p = vcat(p_up, p_down)
    if eltype(z)<:AbstractFloat
        i = unique(i -> z[i], eachindex(z))
    else
        # Fallback for AD tracing...
        i = eachindex(z)
    end
    return Jutul.LinearInterpolant(z[i], p[i])
end

function integrate_phase_density(z_datum, z_end, p0, density_f, phase; n = 1000, g = Jutul.gravity_constant)
    @assert isfinite(p0) "Pressure at contact must be finite"
    z_datum, z_end, p0, g = promote(z_datum, z_end, p0, g)
    T = typeof(p0)
    dz = (z_end - z_datum)/n
    pressure = zeros(T, n+1)
    z = zeros(T, n+1)
    pressure[1] = p0
    z[1] = z_datum

    for i in 2:(n+1)
        p = pressure[i-1]
        depth = z[i-1] + dz
        dens = density_f(p, depth, phase)
        pressure[i] = p + dz*dens*g
        @assert isfinite(pressure[i]) "Equilibriation returned non-finite pressures."
        z[i] = depth
    end
    return (z, pressure)
end

function determine_saturations(depths, contacts, pressures; s_min = missing, s_max = missing, pc = missing)
    nc = length(depths)
    T = promote_type(eltype(depths), eltype(pressures), eltype(contacts))
    nph = length(contacts) + 1
    if ismissing(s_min)
        s_min = [zeros(T, nc) for i in 1:nph]
    end
    if ismissing(s_max)
        s_max = [ones(T, nc) for i in 1:nph]
    end
    sat = zeros(T, nph, nc)
    sat_pc = similar(sat)
    if isnothing(pc) || ismissing(pc)
        for i in eachindex(depths)
            z = depths[i]
            ph = current_phase_index(z, contacts)
            for j in axes(sat, 1)
                is_main = ph == j
                s = is_main*s_max[j][i] + !is_main*s_min[j][i]
                sat[j, i] = s
            end
        end
    else
        ref_ix = 2
        offset = 1
        for ph in 1:nph
            if ph != ref_ix
                s, pc_pair = pc[offset]
                pc_max = maximum(pc_pair)
                pc_min = minimum(pc_pair)
                I = get_1d_interpolator(pc_pair, s, constant_dx = false)
                I_pc = get_1d_interpolator(s, pc_pair, constant_dx = false)
                for i in eachindex(depths)
                    z = depths[i]

                    dp = pressures[ph, i] - pressures[ref_ix, i]
                    if dp > pc_max
                        s_eff = s_max[ph][i]
                    elseif dp < pc_min
                        s_eff = s_min[ph][i]
                    else
                        s_eff = I(dp)
                    end
                    s_eff = clamp(s_eff, s_min[ph][i], s_max[ph][i])
                    sat[ph, i] = s_eff
                    sat_pc[ph, i] = I_pc(s_eff)
                end
                offset += 1
            end
        end
        bad_cells = Int[]
        for i in eachindex(depths)
            sat_i = view(sat, :, i)
            sat_tot = sum(sat_i)
            s_fill = 1 - sat_tot
            if s_fill < 0
                push!(bad_cells, i)
                sat[ref_ix, i] = 0.0
                for j in axes(sat, 1)
                    sat[j, i] /= sat_tot
                end
            else
                sat[ref_ix, i] = s_fill
            end
        end
        if length(bad_cells) > 0
            jutul_message("Initialization", "Negative saturation in $(length(bad_cells)) cells for phase $ref_ix. Normalizing.", color = :yellow)
        end
    end
    return (sat, sat_pc)
end

function current_phase_index(z, depths; reverse = true)
    n = length(depths)+1
    out = -1
    if reverse
        i = current_phase_index(z, Base.reverse(depths), reverse = false)
        out = n - i + 1
    else
        if z < depths[1]
            out = 1
        elseif z > depths[end]
            out = n
        else
            for (i, d) in enumerate(depths)
                if d >= z
                    out = i
                    break
                end
            end
        end
    end
    @assert out > 0
    return out
end
