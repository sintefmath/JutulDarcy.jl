function setup_wells(d::AFIInputFile, reservoir; perf_sort = Dict())
    welldefs = find_records(d, "WellDef", "IX", steps = true, model = true)
    well_dict = OrderedDict{String, Any}()
    for welldef in welldefs
        w2c = get(welldef.value, "WellToCellConnections", Dict())
        if !haskey(w2c, "Cell")
            continue
        end
        wname = welldef.value["WellName"]
        if haskey(well_dict, wname)
            current_w2c = well_dict[wname]["w2c"]
            for (k, v) in pairs(w2c)
                if k in ("Transmissibility", "Status", "PiMultiplier")
                    # These can vary over time - we skip these
                    continue
                end
                if haskey(current_w2c, k)
                    old_value = current_w2c[k]
                    new_value = v
                    if old_value != new_value
                        # Sometimes these get set to zero, which we allow
                        if k == "WellBoreRadius" && length(old_value) == length(new_value)
                            ok = true
                            for i in eachindex(old_value, new_value)
                                ov = old_value[i]
                                nv = new_value[i]
                                if !isapprox(ov, nv, atol = 1e-3) && nv > 0.0
                                    ok = false
                                    break
                                end
                            end
                            if ok
                                continue
                            end
                        end
                        msg = ""
                        for (i, v) in enumerate(new_value)
                            if v isa Real && !isapprox(v, old_value[i], rtol = 1e-6) && !(k == "WellBoreRadius" && v ≈ 0.0)
                                msg *= " Index $i: old=$(old_value[i]) vs new=$v\n"
                            elseif v isa String && v != old_value[i]
                                msg *= " Index $i: old=$(old_value[i]) vs new=$v\n"
                            end
                        end
                        @warn "Duplicate WellDef WellToCellConnections entry for well $wname entry in AFI file with different values. Keyword: $k. Using initially provided entry. Details:\n$msg"
                    end
                else
                    current_w2c[k] = v
                end
            end
        else
            well_dict[wname] = Dict()
            well_dict[wname]["w2c"] = w2c
            well_dict[wname]["ref_depth"] = nothing
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
        cells = get(w2c, "Cell", Int[])
        active = Bool[]
        cells_mapped = Int[]
        head = missing
        for c in cells
            if c isa Int
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
                    head = cell_ijk(mesh, c)
                end
                push!(active, true)
                push!(cells_mapped, next)
            end
        end
        tvd = get(w2c, "TrueVerticalDepth", missing)
        if !ismissing(tvd)
            if isnothing(v["ref_depth"]) && length(cells) > 0
                v["ref_depth"] = w2c["TrueVerticalDepth"][1]
            end
            tvd = tvd[active]
        end
        nperf = length(cells)
        skin = get(w2c, "Skin", zeros(nperf))[active]
        pi_mult = get(w2c, "PiMultiplier", ones(nperf))[active]
        dir = get(w2c, "PenetrationDirection", missing)
        if ismissing(dir)
            dir = fill(:z, nperf)
        else
            dir = map(x -> dir_to_symbol[x], dir)
        end
        dir = dir[active]
        WI = get_if_active(w2c, "Transmissibility", active)
        Kh = get_if_active(w2c, "PermeabilityThickness", active)
        drainage_radius = get_if_active(w2c, "PressureEquivalentRadius", active)
        r = get(w2c, "WellBoreRadius", fill(0.1, nperf))[active]
        compnames = get(w2c, "Completion", missing)
        if ismissing(compnames)
            compnames = map(i -> "COMPLETION_$i", cells)
        end
        compnames = compnames[active]
        worder = get(perf_sort, k, :track)
        sorted_ix = JutulDarcy.well_completion_sortperm(reservoir, head, worder, cells_mapped, dir)
        maybe_sort(::Missing) = missing
        maybe_sort(arr::AbstractVector) = arr[sorted_ix]

        cells_mapped = cells_mapped[sorted_ix]
        skin = maybe_sort(skin)
        dir = maybe_sort(dir)
        WI = maybe_sort(WI)
        Kh = maybe_sort(Kh)
        r = maybe_sort(r)
        tvd = maybe_sort(tvd)
        pi_mult = maybe_sort(pi_mult)

        compnames = maybe_sort(compnames)
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
        w[:well_index_multiplier, JutulDarcy.Perforations()] = Float64.(pi_mult)
        w[:perforation_names, JutulDarcy.Perforations()] = compnames
        w[:original_perforation_indices, JutulDarcy.Perforations()] = sorted_ix
        if !ismissing(tvd)
            w[:perforation_centroids, JutulDarcy.Perforations()][3, :] .= tvd
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
