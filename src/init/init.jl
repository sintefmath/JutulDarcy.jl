
function equilibriate_state(model, contacts, datum_depth = missing, datum_pressure = JutulDarcy.DEFAULT_MINIMUM_PRESSURE;
    cells = missing,
    rs = missing,
    kwarg...)
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
        datum_depth = pmin
    end
    sys = model.system

    init = Dict{Symbol, Any}()
    init = equilibriate_state!(init, pts, model, sys, contacts, datum_depth, datum_pressure; cells = cells, rs = rs, kwarg...)

    is_blackoil = sys isa StandardBlackOilSystem
    if is_blackoil
        sat = init[:Saturations]
        pressure = init[:Pressure]
        nph, nc = size(sat)
        if ismissing(rs)
            Rs = zeros(nc)
        else
            Rs = rs.(pts)
        end
        # TODO: Handle Rv.
        Rv = zeros(nc)
        init[:Rs] = Rs
    end
    return init
end

function equilibriate_state!(init, depths, model, sys, contacts, depth, datum_pressure;
        cells = 1:length(depths),
        rs = missing,
        contacts_pc = missing,
        kwarg...
    )
    if ismissing(contacts_pc)
        contacts_pc = zeros(number_of_phases(sys)-1)
    end
    zmin = minimum(depths)
    zmax = maximum(depths)

    nph = number_of_phases(sys)

    @assert length(contacts) == nph-1

    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)
    disgas = JutulDarcy.has_disgas(sys)
    if disgas
        if JutulDarcy.has_other_phase(sys)
            _, rhoOS, rhoGS = rho_s
        else
            rhoOS, rhoGS = rho_s
        end
    end
    pvt = model.secondary_variables[:PhaseMassDensities].pvt
    relperm = model.secondary_variables[:RelativePermeabilities]
    function density_f(p, z, ph)
        if phases[ph] == LiquidPhase() && disgas
            Rs = min(rs(z), sys.rs_max(p))
            b = JutulDarcy.shrinkage(pvt[ph], 1, p, Rs, 1)
            rho = b*(rhoOS + Rs*rhoGS)
        else
            rho = rho_s[ph]*JutulDarcy.shrinkage(only(pvt[ph].tab), p)
        end
        return rho
    end
    pressures = determine_hydrostatic_pressures(depths, depth, zmin, zmax, contacts, datum_pressure, density_f, contacts_pc)
    s, pc = determine_saturations(depths, contacts, pressures; kwarg...)

    kr = similar(s)
    JutulDarcy.update_kr!(kr, relperm, model, s, cells)

    init[:Saturations] = s
    init[:Pressure] = init_reference_pressure(pressures, contacts, pc, kr, 2)
    return init
end


function parse_state0_equil(model, datafile)
    sys = model.system
    init = Dict{Symbol, Any}()
    sol = datafile["SOLUTION"]
    G = physical_representation(model.data_domain)
    nc = number_of_cells(G)
    nph = number_of_phases(model.system)
    ix = G.cell_map
    is_blackoil = sys isa StandardBlackOilSystem
    has_water = haskey(datafile["RUNSPEC"], "WATER")
    has_oil = haskey(datafile["RUNSPEC"], "OIL")
    has_gas = haskey(datafile["RUNSPEC"], "GAS")
    disgas = JutulDarcy.has_disgas(model.system)

    equil = sol["EQUIL"]
    nequil = JutulDarcy.InputParser.number_of_tables(datafile, :equil)
    @assert length(equil) == nequil
    inits = []
    for ereg in 1:nequil
        eq = equil[ereg]
        cells = findall(isequal(ereg), model.data_domain[:eqlnum])
        datum_depth = eq[1]
        datum_pressure = eq[2]

        woc = eq[3]
        woc_pc = eq[4]
        goc = eq[5]
        goc_pc = eq[6]
        # Contact depths
        s_max = 1.0
        s_min = 0.0

        non_connate = 1.0
        s_max = Float64[]
        s_min = Float64[]
        # @assert !haskey(model.secondary_variables, :CapillaryPressure) "Capillary initialization not yet implemented."
        has_pc = haskey(model.secondary_variables, :CapillaryPressure)
        if has_pc
            pc_f = model.secondary_variables[:CapillaryPressure].pc
            pc = []
            for (i, f) in enumerate(pc_f)
                s = only(f).X
                cap = only(f).F
                if s[1] < 0
                    s = s[2:end]
                    cap = cap[2:end]
                end
                ix = unique(i -> cap[i], 1:length(cap))

                if nph == 3 && i == 1
                    @. cap *= -1
                end
                push!(pc, (s = s[ix], pc = cap[ix]))
            end
        else
            pc = nothing
        end
        if has_water
            krw = only(model.secondary_variables[:RelativePermeabilities].krw)
            swcon = krw.critical
            push!(s_min, swcon)
            push!(s_max, non_connate)
            non_connate = 1.0 - swcon
        end
        if has_oil
            push!(s_min, 0.0)
            push!(s_max, non_connate)
        end
        if has_gas
            push!(s_min, 0.0)
            push!(s_max, non_connate)
        end

        if nph == 1
            error("Not implemented.")
        elseif nph == 2
            if has_oil && has_gas
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

        if disgas
            @assert haskey(sol, "RSVD")
            rsvd = sol["RSVD"][ereg]
            z = rsvd[:, 1]
            Rs = rsvd[:, 2]
            rs = Jutul.LinearInterpolant(z, Rs)
        else
            rs = missing
        end

        subinit = equilibriate_state(model, contacts, datum_depth, datum_pressure, cells = cells, contacts_pc = contacts_pc, s_min = s_min, s_max = s_max, rs = rs, pc = pc)
        push!(inits, subinit)
    end
    # TODO: Handle multiple regions by merging each init
    return only(inits)
end

function init_reference_pressure(pressures, contacts, kr, pc, ref_ix = 2)
    nph, nc = size(kr)
    p = zeros(nc)
    for i in eachindex(p)
        p[i] = pressures[ref_ix, i]
        if kr[ref_ix, i] <= 0
            for ph in 1:nph
                if kr[ph, i] > 0
                    p[i] = pressures[ph, i]# - pc[ph, i]
                end
            end
        end
    end
    return p
end

function determine_hydrostatic_pressures(depths, depth, zmin, zmax, contacts, datum_pressure, density_f, contacts_pc)
    nc = length(depths)
    nph = length(contacts) + 1
    ref_ix = 2
    I_ref = phase_pressure_depth_table(depth, zmin, zmax, datum_pressure, density_f, ref_ix)
    pressures = zeros(nph, nc)
    pos = 1
    for ph in 1:nph
        if ph == ref_ix
            I = I_ref
        else
            contact = contacts[pos]
            datum_pressure_ph = I_ref(contact)
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
    return Jutul.LinearInterpolant(z, p)
end

function integrate_phase_density(z_datum, z_end, p0, density_f, phase; n = 1000, g = Jutul.gravity_constant)
    dz = (z_end - z_datum)/n
    pressure = zeros(n+1)
    z = zeros(n+1)
    pressure[1] = p0
    z[1] = z_datum

    for i in 2:(n+1)
        p = pressure[i-1]
        depth = z[i-1] + dz
        pressure[i] = p + dz*density_f(p, depth, phase)*g
        z[i] = depth
    end
    return (z, pressure)
end

function determine_saturations(depths, contacts, pressures; s_min = missing, s_max = missing, pc = nothing)
    nc = length(depths)
    nph = length(contacts) + 1
    if ismissing(s_max)
        s_max = zeros(nph)
    end
    if ismissing(s_max)
        s_max = zeros(nph)
    end
    sat = zeros(nph, nc)
    sat_pc = similar(sat)
    if isnothing(pc)
        for i in eachindex(depths)
            z = depths[i]
            ph = current_phase_index(z, contacts)
            for j in axes(sat, 1)
                is_main = ph == j
                s = is_main*s_max[j] + !is_main*s_min[j]
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

                I = Jutul.LinearInterpolant(pc_pair, s)
                I_pc = Jutul.LinearInterpolant(s, pc_pair)
                for i in eachindex(depths)
                    z = depths[i]

                    dp = pressures[ph, i] - pressures[ref_ix, i]
                    if dp > pc_max
                        s_eff = s_max[ph]
                    elseif dp < pc_min
                        s_eff = s_min[ph]
                    else
                        s_eff = I(dp)
                    end
                    s_eff = clamp(s_eff, s_min[ph], s_max[ph])
                    sat[ph, i] = s_eff
                    sat_pc[ph, i] = I_pc(s_eff)
                end
                offset += 1
            end
        end
        for i in eachindex(depths)
            s_fill = 1 - sum(view(sat, :, i))
            if s_fill < 0
                @warn "Negative saturation in cell $i: $s_fill"
            end
            sat[ref_ix, i] = s_fill
        end
    end
    return (sat, sat_pc)
end

function current_phase_index(z, depths; reverse = true)
    n = length(depths)+1
    if reverse
        i = current_phase_index(z, Base.reverse(depths), reverse = false)
        return n - i + 1
    else
        if z < depths[1]
            return 1
        elseif z > depths[end]
            return n
        else
            for (i, d) in enumerate(depths)
                if d > z
                    return i
                end
            end
        end
    end
end
