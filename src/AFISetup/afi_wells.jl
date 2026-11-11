function afi_completion_names(welldefs)
    names = Dict{Tuple{String, Any}, Union{String, Missing}}()
    for record in welldefs
        wname = record.value["WellName"]
        connections = get(record.value, "WellToCellConnections", Dict())
        cells = get(connections, "Cell", missing)
        completions = get(connections, "Completion", missing)
        if ismissing(cells) || ismissing(completions) || isempty(completions)
            continue
        end
        length(completions) == length(cells) ||
            throw(ArgumentError("Well $wname has different numbers of cells and completion names"))
        for (cell, completion) in zip(cells, completions)
            ismissing(completion) && continue
            key = (wname, cell)
            name = String(completion)
            if !isequal(get!(names, key, name), name)
                names[key] = missing
            end
        end
    end
    return names
end

function afi_connection_ids(wname, connections, names)
    cells = get(connections, "Cell", missing)
    ismissing(cells) && return nothing
    completions = get(connections, "Completion", missing)
    ids = map(eachindex(cells)) do i
        cell = cells[i]
        if ismissing(completions) || isempty(completions) || ismissing(completions[i])
            name = get(names, (wname, cell), nothing)
            ismissing(name) &&
                throw(ArgumentError("Well $wname has multiple completions in cell $cell. The AFI record must identify the completion."))
        else
            name = String(completions[i])
        end
        return (cell, name)
    end
    allunique(ids) || throw(ArgumentError("Well $wname has duplicate cell and completion pairs in a WellDef record"))
    return ids
end

function merge_afi_perforation_properties!(properties, connections, positions, wname)
    for field in ("TrueVerticalDepth", "Skin", "PiMultiplier", "PenetrationDirection",
            "Transmissibility", "PermeabilityThickness", "WellBoreRadius")
        values = get(connections, field, missing)
        if ismissing(values) || isempty(values)
            continue
        end
        if field != "PiMultiplier" && length(values) != length(positions)
            throw(ArgumentError("Well $wname has different numbers of cells and $field values"))
        end
        column = get!(Dict{Int, Any}, properties, field)
        changes_over_time = field in ("PiMultiplier", "Transmissibility")
        for (row, value) in enumerate(values)
            ismissing(value) && continue
            row <= length(positions) ||
                throw(ArgumentError("Well $wname has a $field value at row $row with no corresponding AFI connection"))
            position = positions[row]
            old = get(column, position, missing)
            if ismissing(old) || (!changes_over_time && old isa Real && old ≈ 0.0)
                column[position] = value
            elseif !changes_over_time
                if field == "WellBoreRadius"
                    matches = value ≈ 0.0 || isapprox(old, value; atol = 1e-3)
                elseif old isa Real && value isa Real
                    matches = isapprox(old, value; rtol = 1e-6)
                else
                    matches = isequal(old, value)
                end
                if !matches
                    @warn "Conflicting $field for well $wname, connection $(connections["Cell"][row]). Using the first value ($old)."
                end
            end
        end
    end
end

function setup_wells(d::AFIInputFile, reservoir; perf_sort = :input)
    welldefs = find_records(d, "WellDef", "IX", steps = true, model = true)
    completion_names = afi_completion_names(welldefs)
    well_dict = OrderedDict{String, Any}()
    for record in welldefs
        wname = record.value["WellName"]
        connections = get(record.value, "WellToCellConnections", Dict())
        well = get!(well_dict, wname) do
            Dict{String, Any}(
                "positions" => OrderedDict{Tuple, Int}(),
                "properties" => Dict{String, Dict{Int, Any}}(),
                "ref_depth" => nothing
            )
        end
        ids = afi_connection_ids(wname, connections, completion_names)
        (isnothing(ids) || isempty(ids)) && continue
        positions = well["positions"]
        rows = [get!(() -> length(positions) + 1, positions, id) for id in ids]
        merge_afi_perforation_properties!(well["properties"], connections, rows, wname)
    end
    # We do a second pass when all wells have been found
    well_kws = find_records(d, "Well", "FM", steps = true, model = true)
    for well_kw in well_kws
        val = well_kw.value
        if haskey(val, "BottomHoleRefDepth")
            welltab = Dict(well_kw.value["name"] => val)
        else
            welltab = get(val, "wells", Dict())
        end
        for (wname, val) in pairs(welltab)
            if haskey(val, "BottomHoleRefDepth") && haskey(well_dict, wname)
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
    dir_to_symbol = Dict(
        GeoEnergyIO.IXParser.IX_I => :x,
        GeoEnergyIO.IXParser.IX_J => :y,
        GeoEnergyIO.IXParser.IX_K => :z,
    )
    mesh = physical_representation(reservoir)
    ijk_to_local = Dict(cell_ijk(mesh, i) => i for i in 1:number_of_cells(mesh))
    global_to_local = afi_to_jutul_cell_map(d, mesh)
    wells = []
    for (name, well) in pairs(well_dict)
        ids = collect(keys(well["positions"]))
        isempty(ids) && continue
        properties = well["properties"]
        input_indices = Int[]
        local_cells = Int[]
        for (i, (cell, _)) in enumerate(ids)
            local_cell = get(cell isa Int ? global_to_local : ijk_to_local, cell, nothing)
            if !isnothing(local_cell)
                push!(input_indices, i)
                push!(local_cells, local_cell)
            end
        end
        if isempty(local_cells)
            @warn "Skipping AFI well $name because it has no active connections"
            continue
        end

        dirs = get(properties, "PenetrationDirection", Dict())
        direction = Symbol[dir_to_symbol[get(dirs, i, GeoEnergyIO.IXParser.IX_K)] for i in input_indices]
        head = cell_ijk(mesh, first(local_cells))
        order = perf_sort isa Symbol ? perf_sort : get(perf_sort, name, :input)
        permutation = JutulDarcy.well_completion_sortperm(reservoir, head, order, local_cells, direction)
        input_indices = input_indices[permutation]
        local_cells = local_cells[permutation]

        function property_values(field, default)
            column = get(properties, field, Dict())
            return [get(column, i, default) for i in input_indices]
        end
        depth = get(properties, "TrueVerticalDepth", Dict())
        reference_depth = well["ref_depth"]
        if isnothing(reference_depth)
            reference_depth = get(depth, 1, nothing)
        end
        w = setup_well(reservoir, local_cells,
            skin = Float64.(property_values("Skin", 0.0)),
            dir = direction[permutation],
            WI = Float64.(property_values("Transmissibility", NaN)),
            Kh = Float64.(property_values("PermeabilityThickness", NaN)),
            radius = Float64.(property_values("WellBoreRadius", 0.1)),
            reference_depth = reference_depth,
            name = Symbol(name)
        )
        perfs = JutulDarcy.Perforations()
        connection_ids = ids[input_indices]
        w[:well_index_multiplier, perfs] = Float64.(property_values("PiMultiplier", 1.0))
        w[:perforation_names, perfs] = [isnothing(completion) ? "COMPLETION_$cell" : completion for (cell, completion) in connection_ids]
        w[:original_perforation_indices, perfs] = input_indices
        w[:perforation_connection_ids, perfs] = connection_ids
        if !isempty(depth)
            centroids = w[:perforation_centroids, perfs]
            for (i, input_index) in enumerate(input_indices)
                centroids[3, i] = get(depth, input_index, centroids[3, i])
            end
        end
        push!(wells, w)
    end
    return wells
end

function afi_to_jutul_cell_map(d, mesh)
    cell_offset = cell_index_offset(d)
    global_to_local = Dict{Int, Int}()
    cmap = mesh isa UnstructuredMesh ? mesh.cell_map : nothing
    if isnothing(cmap)
        cmap = 1:number_of_cells(mesh)
    end
    for (i, c) in enumerate(cmap)
        c = c - 1 + cell_offset
        global_to_local[c] = i
    end
    return global_to_local
end
