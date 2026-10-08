function gas_water_component_index(eos)
    w = findfirst(isequal(MultiComponentFlash.COMPONENT_H2O), eos.type.component_types)
    isnothing(w) && throw(ArgumentError("GASWAT requires a water component."))
    return w
end

# EQUIL mode 2 supplies a vapor composition at the contact. At other pressures,
# adjust the trial feed between ordinary flashes to preserve its dry-gas ratios.
function gas_water_phase_equilibrium(eos, p, T, gas_composition)
    MCF = MultiComponentFlash
    n = MCF.number_of_components(eos)
    w = gas_water_component_index(eos)
    y = Vector{Float64}(gas_composition)
    y /= sum(y)
    dry = copy(y)
    dry[w] = 0.0
    dry /= sum(dry)
    x = zeros(n)
    x[w] = 1.0
    for it in 1:100
        # Equal phase mole amounts keep the trial feed inside the two-phase region.
        z = (x + y)/2
        # Keep absent components below the composition-constraint tolerance.
        f = MCF.flashed_mixture_2ph(eos, (p = p, T = T, z = z);
            tolerance = 1e-10, maxiter = 200, z_min = eps(Float64))
        0 < f.V < 1 || error("Expected two-phase gas-water equilibrium at p=$p, T=$T.")
        x = f.liquid.mole_fractions
        y = f.vapor.mole_fractions
        ynext = (1.0 - y[w])*dry
        ynext[w] = y[w]
        if maximum(abs, ynext - y) < 1e-12
            return f
        end
        y = ynext
    end
    error("Gas-water composition constraint failed to converge at p=$p, T=$T.")
end

function gas_water_dew_pressure(eos, T, y)
    w = gas_water_component_index(eos)
    yw = y[w]/sum(y)
    residual(p) = gas_water_phase_equilibrium(eos, p, T, y).vapor.mole_fractions[w] - yw
    lo = 1e5
    flo = residual(lo)
    hi = lo
    for it in 1:40
        hi *= 1.4
        fhi = residual(hi)
        if flo*fhi <= 0.0
            for _ in 1:70
                mid = (lo + hi)/2
                fm = residual(mid)
                if abs(hi - lo) < 1e-8*mid
                    return mid
                elseif flo*fm <= 0.0
                    hi = mid
                else
                    lo, flo = mid, fm
                end
            end
        end
        lo, flo = hi, fhi
        hi < 1e8 || break
    end
    throw(ArgumentError("No GASWAT dew pressure found for the specified vapor composition at T=$T."))
end

function gas_water_composition_function(table, n)
    depths = table[:, 1]
    values = [SVector{n}(table[i, 2:n+1]) for i in axes(table, 1)]
    values = map(x -> x/sum(x), values)
    return get_1d_interpolator(depths, values)
end

function gas_water_capillary_table(model, props, sreg)
    if !haskey(model.secondary_variables, :CapillaryPressure)
        return nothing
    end
    table = props["GSF"][sreg]
    indices = unique(i -> table[i, 3], axes(table, 1))
    saturation, capillary = table[indices, 1], table[indices, 3]
    # Keep the terminal plateau and its maximum saturation when inverting the
    # table. One ulp makes the inverse well-defined without changing the ramp.
    if saturation[end] < table[end, 1]
        push!(saturation, table[end, 1])
        push!(capillary, nextfloat(capillary[end]))
    end
    return [(s = saturation, pc = capillary)]
end

function gas_water_overall_composition(eos, p, T, saturation, equilibrium)
    MCF = MultiComponentFlash
    liquid = equilibrium.liquid
    vapor = equilibrium.vapor
    nl = saturation[1]/MCF.molar_volume(eos, p, T, liquid)
    nv = saturation[2]/MCF.molar_volume(eos, p, T, vapor)
    return (nl*liquid.mole_fractions + nv*vapor.mole_fractions)/(nl + nv)
end

function parse_state0_gas_water_equil(model, datafile; cell_nz = 1)
    MCF = MultiComponentFlash
    sys = model.system
    eos = sys.equation_of_state
    domain = model.data_domain
    props, sol = datafile["PROPS"], datafile["SOLUTION"]
    n = MCF.number_of_components(eos)
    nc = number_of_cells(domain)
    eqlnum = domain[:eqlnum]
    satnum = domain[:satnum]
    init = Dict{Symbol, Any}(:Pressure => zeros(nc),
        :OverallMoleFractions => zeros(n, nc), :Temperature => zeros(nc),
        :Saturations => zeros(2, nc))
    depths = domain[:cell_centroids][3, :]
    kr = model.secondary_variables[:RelativePermeabilities]
    for (ereg, eq) in enumerate(sol["EQUIL"])
        eq[10] == 2 || throw(ArgumentError("GASWAT EQUIL currently requires vapor composition at the contact (item 10 = 2)."))
        datum_depth, datum_pressure, contact, contact_pc = eq[1:4]
        datum_depth == contact || throw(ArgumentError("GASWAT EQUIL item 10 = 2 requires the datum at the gas-water contact."))
        composition = gas_water_composition_function(props["ZMFVD"][ereg], n)
        T_z = equil_temperature_function(datafile, ereg)
        if ismissing(T_z)
            throw(ArgumentError("GASWAT EQUIL needs TEMPVD or RTEMP."))
        end
        if eq[11] != 1 || !isfinite(datum_pressure)
            datum_pressure = gas_water_dew_pressure(eos, T_z(contact), composition(contact))
        end
        # The specified dew pressure belongs to the vapor; liquid is the model's
        # reference phase. Capillary pressure is p_gas - p_water.
        datum_pressure -= contact_pc
        function density_function(p, z, phase)
            f = gas_water_phase_equilibrium(eos, p, T_z(z), composition(z))
            if phase == 1
                ph = MCF.phase_data(f, Val(:liquid))
            else
                ph = MCF.phase_data(f, Val(:vapor))
            end
            return MCF.mass_density(eos, p, T_z(z), ph)
        end
        cells_equil = findall(isequal(ereg), eqlnum)
        for sreg in unique(satnum[cells_equil])
            cells_sat = findall(isequal(sreg), satnum)
            cells = intersect_sorted(cells_equil, cells_sat)
            krw = table_by_region(kr.krw, sreg)
            krg = table_by_region(kr.krg, sreg)
            smin = [fill(krw.connate, length(cells)), fill(krg.connate, length(cells))]
            smax = [fill(krw.input_s_max, length(cells)), fill(krg.input_s_max, length(cells))]
            pc = gas_water_capillary_table(model, props, sreg)
            subinit = equilibriate_state(model, (contact,), datum_depth, datum_pressure;
                cells = cells, contacts_pc = (contact_pc,),
                composition = composition, T_z = T_z,
                density_function = density_function, s_min = smin, s_max = smax,
                satnum = sreg, cell_nz = cell_nz, pc = pc)
            fill_subinit!(init[:Pressure], cells, subinit[:Pressure])
            fill_subinit!(init[:Saturations], cells, subinit[:Saturations])
            for (i, c) in enumerate(cells)
                p = subinit[:Pressure][i]
                T = T_z(depths[c])
                f = gas_water_phase_equilibrium(eos, p, T, composition(depths[c]))
                # Convert the hydrostatic phase volumes to component mole
                # fractions. Flashing this mixture must recover the saturations.
                saturation = view(subinit[:Saturations], :, i)
                z = gas_water_overall_composition(eos, p, T, saturation, f)
                init[:Temperature][c] = T
                init[:OverallMoleFractions][:, c] .= z
            end
        end
    end
    return init
end
