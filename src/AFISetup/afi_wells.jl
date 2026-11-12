# Use a name supplied elsewhere only when the cell has a single named completion.
function afi_completion_names(welldefs)
    names = Dict{Tuple{String, Any}, Union{String, Missing}}()
    for record in welldefs
        w2c = get(record.value, "WellToCellConnections", Dict())
        cells = get(w2c, "Cell", missing)
        completion = get(w2c, "Completion", missing)
        (ismissing(cells) || ismissing(completion)) && continue
        for (cell, name) in zip(cells, completion)
            ismissing(name) && continue
            key = (record.value["WellName"], cell)
            name = String(name)
            isequal(get!(names, key, name), name) || (names[key] = missing)
        end
    end
    return names
end

function afi_connection_ids(wname, w2c, names)
    cells = get(w2c, "Cell", missing)
    ismissing(cells) && return nothing
    completion = get(w2c, "Completion", missing)
    return map(eachindex(cells)) do i
        cell = cells[i]
        name = ismissing(completion) ? missing : get(completion, i, missing)
        return (cell, ismissing(name) ? get(names, (wname, cell), nothing) : String(name))
    end
end

function setup_wells(d::AFIInputFile, reservoir; perf_sort = Dict())
    welldefs = find_records(d, "WellDef", "IX", steps = true, model = true)
    names = afi_completion_names(welldefs)
    well_dict = OrderedDict{String, Any}()
    for welldef in welldefs
        w2c = get(welldef.value, "WellToCellConnections", Dict())
        wname = welldef.value["WellName"]
        ids = afi_connection_ids(wname, w2c, names)
        (isnothing(ids) || isempty(ids)) && continue
        first_record = !haskey(well_dict, wname)
        well = get!(well_dict, wname) do
            Dict{String, Any}("w2c" => Dict{String, Any}("UniqueID" => Tuple[]), "ref_depth" => nothing)
        end
        current_w2c = well["w2c"]
        current_uids = current_w2c["UniqueID"]
        uid_pos = map(ids) do id
            pos = findfirst(isequal(id), current_uids)
            if isnothing(pos)
                push!(current_uids, id)
                pos = length(current_uids)
            end
            return pos
        end

        for k in union(keys(current_w2c), keys(w2c))
            k in ("UniqueID", "Cell", "Completion", "Status") && continue
            old = get(current_w2c, k, [])
            new = get(w2c, k, [])
            # Keep the initial values for fields that can change in the schedule.
            time_varying = k in ("PiMultiplier", "Transmissibility")
            if time_varying && !first_record
                new = []
            end
            old = ismissing(old) ? [] : old
            new = ismissing(new) ? [] : new
            isempty(old) && isempty(new) && continue
            aligned = Any[get(old, i, missing) for i in eachindex(current_uids)]
            for (i, value) in enumerate(new)
                (i > length(uid_pos) || ismissing(value)) && continue
                pos = uid_pos[i]
                previous = aligned[pos]
                if ismissing(previous) || (!time_varying && previous isa Real && previous ≈ 0.0)
                    aligned[pos] = value
                elseif !time_varying
                    if previous isa Real && value isa Real
                        matches = k == "WellBoreRadius" ?
                            (value ≈ 0.0 || isapprox(previous, value; atol = 1e-3)) :
                            isapprox(previous, value; rtol = 1e-6)
                    else
                        matches = isequal(previous, value)
                    end
                    if !matches
                        @warn "Conflicting $k for well $wname, connection $(ids[i]). Using the first value ($previous)."
                    end
                end
            end
            current_w2c[k] = aligned
        end
    end
    if perf_sort isa Symbol
        perf_sort = Dict{String, Symbol}(
            k => perf_sort for k in keys(well_dict)
        )
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
    dir_to_symbol = Dict(
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
        ids = w2c["UniqueID"]
        cells = first.(ids)
        active = Bool[]
        cells_mapped = Int[]
        head = missing
        for (c, completion) in ids
            if ismissing(completion)
                # An unnamed row cannot identify one of several completions in this cell.
                next = missing
            elseif c isa Int
                # TODO: Not clear to me if these include the cell index offset or not
                # next = c
                next = get(global_to_local, c, missing)
            else
                c::Tuple{Int, Int, Int}
                next = get(ijk_to_local, Tuple(c), missing)
            end
            if ismissing(next)
                push!(active, false)
            else
                if ismissing(head)
                    head = cell_ijk(mesh, next)
                end
                push!(active, true)
                push!(cells_mapped, next)
            end
        end
        isempty(cells_mapped) && continue
        tvd = get_if_active(w2c, "TrueVerticalDepth", active)
        if !ismissing(tvd)
            first_depth = w2c["TrueVerticalDepth"][1]
            if isnothing(v["ref_depth"]) && !ismissing(first_depth)
                v["ref_depth"] = first_depth
            end
        end
        nperf = length(cells)
        skin = get(w2c, "Skin", zeros(nperf))[active]
        pi_mult = get(w2c, "PiMultiplier", ones(nperf))[active]
        dir = get(w2c, "PenetrationDirection", missing)
        if ismissing(dir)
            dir = fill(:z, nperf)
        else
            dir = map(x -> dir_to_symbol[coalesce(x, GeoEnergyIO.IXParser.IX_K)], dir)
        end
        dir = dir[active]
        WI = get_if_active(w2c, "Transmissibility", active)
        Kh = get_if_active(w2c, "PermeabilityThickness", active)
        drainage_radius = get_if_active(w2c, "PressureEquivalentRadius", active)
        r = get(w2c, "WellBoreRadius", fill(0.1, nperf))[active]
        compnames = map(ids[active]) do (cell, name)
            return isnothing(name) ? "COMPLETION_$cell" : name
        end
        worder = get(perf_sort, k, :depth)
        if lowercase(String(worder)) == "depth" && !ismissing(tvd)
            depths = coalesce.(tvd, reservoir[:cell_centroids][3, cells_mapped])
            sorted_ix = sortperm(depths)
        else
            sorted_ix = JutulDarcy.well_completion_sortperm(reservoir, head, worder, cells_mapped, dir)
        end
        function reorder_and_typecheck(::Missing, T = missing, default = missing)
            return missing
        end
        function reorder_and_typecheck(arr::AbstractVector, T = missing, default = missing)
            out = arr[sorted_ix]
            for (i, v) in enumerate(out)
                if ismissing(v)
                    out[i] = default
                end
            end
            if ismissing(T)
                out = out
            else
                out = convert(Vector{T}, out)
            end
            return out
        end

        cells_mapped = cells_mapped[sorted_ix]
        skin = reorder_and_typecheck(skin, Float64, 0.0)
        WI = reorder_and_typecheck(WI, Float64, NaN)
        Kh = reorder_and_typecheck(Kh, Float64, NaN)
        r = reorder_and_typecheck(r, Float64, 0.1)
        pi_mult = reorder_and_typecheck(pi_mult, Float64, 1.0)
        # Not necessarily Float64
        dir = reorder_and_typecheck(dir)
        tvd = reorder_and_typecheck(tvd)
        compnames = reorder_and_typecheck(compnames)
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
        w[:well_index_multiplier, JutulDarcy.Perforations()] = pi_mult
        w[:perforation_names, JutulDarcy.Perforations()] = compnames
        w[:original_perforation_indices, JutulDarcy.Perforations()] = findall(active)[sorted_ix]
        w[:perforation_connection_ids, JutulDarcy.Perforations()] = ids[active][sorted_ix]
        if !ismissing(tvd)
            depths = view(w[:perforation_centroids, JutulDarcy.Perforations()], 3, :)
            depths .= coalesce.(tvd, depths)
        end
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
