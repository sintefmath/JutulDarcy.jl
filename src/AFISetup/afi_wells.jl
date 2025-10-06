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
    end
    dir_to_str = Dict(
        GeoEnergyIO.IXParser.IX_I => :x,
        GeoEnergyIO.IXParser.IX_J => :y,
        GeoEnergyIO.IXParser.IX_K => :z,
    )
    # Set up mappings
    mesh = physical_representation(reservoir)
    # IJK indices
    ijk = map(i -> cell_ijk(mesh, i), 1:number_of_cells(mesh))
    ijk_to_local = Dict{Tuple{Int, Int, Int}, Int}()
    for (i, c) in enumerate(ijk)
        ijk_to_local[c] = i
    end
    # Linear indices
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
    wells = []
    for (k, v) in pairs(well_dict)
        w2c = v["w2c"]
        cells = get(w2c, "Cell", Int[])
        active = Bool[]
        cells_mapped = Int[]
        for c in cells
            if c isa Int
                next = get(global_to_local, c, missing)
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
        WI = get(w2c, "Transmissibility", missing)
        if length(WI) == 0
            WI = missing
        else
            WI = WI[active]
        end
        r = get(w2c, "WellboreRadius", fill(0.1, nperf))[active]
        # TODO: Permeability thickness, multipliers, etc
        w = setup_well(reservoir, cells_mapped, skin = skin, dir = dir, WI = WI, radius = r, name = Symbol(k))
        compnames = get(w2c, "Completion", missing)
        if ismissing(compnames)
            compnames = map(i -> "COMPLETION_$i", cells)
        end
        pprm = w.perforation_parameters
        pprm[:names] = compnames[active]
        push!(wells, w)
    end
    return wells
end
