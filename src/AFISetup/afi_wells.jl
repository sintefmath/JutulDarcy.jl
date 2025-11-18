function setup_wells(d::AFIInputFile, reservoir)
    welldefs = find_records(d, "WellDef", "IX", steps = true, model = true)
    well_dict = OrderedDict{String, Any}()
    for welldef in welldefs
        w2c = get(welldef.value, "WellToCellConnections", Dict())
        if !haskey(w2c, "Cell")
            continue
        end
        wname = welldef.value["WellName"]
        !haskey(well_dict, wname) || error("Duplicate WellDef keyword for well $wname with WellToCellConnections entry in AFI file.")
        well_dict[wname] = Dict()
        well_dict[wname]["w2c"] = w2c
        well_dict[wname]["ref_depth"] = nothing
    end
    # We do a second pass when all wells have been found
    well_kws = find_records(d, "Well", "FM", steps = true, model = true)
    for well_kw in well_kws
        val = well_kw.value
        welltab = get(val, "wells", Dict())
        for (wname, val) in pairs(welltab)
            if haskey(val, "BottomHoleRefDepth")
                newval = val["BottomHoleRefDepth"]
                oldval = well_dict[wname]["ref_depth"]
                if isnothing(oldval)
                    well_dict[wname]["ref_depth"] = newval
                elseif !(oldval ≈ newval)
                    println("Inconsistent BottomHoleRefDepth for well $wname: $oldval vs $newval. Using the first provided value ($oldval).")
                end
            end
        end
    end
    dir_to_str = Dict(
        GeoEnergyIO.IXParser.IX_I => :x,
        GeoEnergyIO.IXParser.IX_J => :y,
        GeoEnergyIO.IXParser.IX_K => :z,
    )
    if haskey(reservoir, :net_to_gross)
        ntg = reservoir[:net_to_gross]
    else
        ntg = missing
    end
    # Set up mappings
    mesh = physical_representation(reservoir)
    # IJK indices
    ijk = map(i -> cell_ijk(mesh, i), 1:number_of_cells(mesh))
    ijk_to_local = Dict{Tuple{Int, Int, Int}, Int}()
    for (i, c) in enumerate(ijk)
        ijk_to_local[c] = i
    end
    # Linear indices
    global_to_local = afi_to_jutul_cell_map(d, mesh)
    function get_if_active(w2c, kw, active)
        x = get(w2c, kw, missing)
        if ismissing(x) || length(x) == 0
            return missing
        end
        return x[active]
    end
    wells = []
    for (k, v) in pairs(well_dict)
        w2c = v["w2c"]
        cells = get(w2c, "Cell", Int[])
        active = Bool[]
        cells_mapped = Int[]
        for c in cells
            if c isa Int
                # TODO: Not clear to me if these include the cell index offset or not
                next = c
                # next = get(global_to_local, c, missing)
            else
                c::Tuple{Int, Int, Int}
                next = get(ijk_to_local, Tuple(c), missing)
            end
            if ismissing(next)
                push!(active, false)
            else
                push!(active, true)
                push!(cells_mapped, next)
            end
        end
        nperf = length(cells)
        skin = get(w2c, "SkinFactors", zeros(nperf))[active]
        dir = get(w2c, "PenetrationDirection", missing)
        if ismissing(dir)
            dir = fill(:z, nperf)
        else
            dir = map(x -> dir_to_str[x], dir)
        end
        dir = dir[active]
        WI = get_if_active(w2c, "Transmissibility", active)
        Kh = get_if_active(w2c, "PermeabilityThickness", active)
        drainage_radius = get_if_active(w2c, "PressureEquivalentRadius", active)
        r = get(w2c, "WellBoreRadius", fill(0.1, nperf))[active]
        w = setup_well(reservoir, cells_mapped,
            skin = skin,
            dir = dir,
            WI = WI,
            Kh = Kh,
            net_to_gross = ntg,
            # drainage_radius = drainage_radius,
            radius = r,
            reference_depth = v["ref_depth"],
            name = Symbol(k)
        )
        compnames = get(w2c, "Completion", missing)
        if ismissing(compnames)
            compnames = map(i -> "COMPLETION_$i", cells)
        end
        w[:perforation_names, JutulDarcy.Perforations()] = compnames[active]
        push!(wells, w)
    end
    return wells
end

function afi_to_jutul_cell_map(d, mesh)
    cell_offset = cell_index_offset(d)
    global_to_local = Dict{Int, Int}()
    cmap = mesh.cell_map
    if isnothing(cmap)
        cmap = 1:number_of_cells(mesh)
    end
    for (i, c) in enumerate(cmap)
        c = c - 1 + cell_offset
        global_to_local[c] = i
    end
    return global_to_local
end
