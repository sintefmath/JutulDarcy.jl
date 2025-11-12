function JutulDarcy.reservoir_domain(d::AFIInputFile; mesh = missing, kwarg...)
    mesh = setup_mesh_afi(d, mesh)

    # TODO: Handle regions and better warnings
    return setup_reservoir_domain_afi(d, mesh; kwarg...)
end

function setup_mesh_afi(d::AFIInputFile, mesh)
    return mesh
end

function cartdims_from_structured_info(d::AFIInputFile)
    sinfo = find_records(d, "StructuredInfo", once = true)
    !isnothing(sinfo) || error("No StructuredInfo record found in AFI file.")
    d = sinfo.value
    return (d["I"], d["J"], d["K"])
end

function cell_index_offset(d::AFIInputFile)
    sinfo = find_records(d, "StructuredInfo", once = true)
    start = 0
    if !isnothing(sinfo)
        d = sinfo.value
        start = get(d, "FirstCellId", start)
    end
    return start
end

function setup_mesh_afi(afi::AFIInputFile)
    return setup_mesh_afi(afi, missing)
end

function setup_mesh_afi(afi::AFIInputFile, mesh::Missing)
    IX = afi.setup["IX"]
    if haskey(IX, "RESQML")
        resqml = IX["RESQML"]
        haskey(IX["RESQML"], "GRID") || error("No GRID section in converted AFI file under the RESQML section.")
        mesh = GeoEnergyIO.mesh_from_grid_section(resqml["GRID"])
    else
        pillar_grid = find_records(afi, "StraightPillarGrid", once = true)
        if isnothing(pillar_grid)
            error("No RESQML section in AFI file, and no StraightPillarGrid record found.")
        end
        pg = pillar_grid.value
        cartdims = cartdims_from_structured_info(afi)
        nx, ny, nz = cartdims
        dx = pg["DeltaX"]
        dy = pg["DeltaY"]
        dz = pg["DeltaZ"]
        tops = pg["PillarTops"]

        dcp = pg["CellDoubleProperty"]
        if ismissing(tops) && haskey(dcp, "CELL_TOP_DEPTH")
            tops = dcp["CELL_TOP_DEPTH"]
        end
        if !ismissing(tops)
            tops = tops[1:nx*ny]
        end
        grid = Dict(
            "cartDims" => cartdims,
            "DXV" => dx,
            "DYV" => dy,
            "DZV" => dz,
            "TOPS" => tops,
        )
        mesh = mesh_from_grid_section(grid)
    end
    !ismissing(mesh) || error("Could not set up mesh from AFI file.")
    return mesh
end

function setup_reservoir_domain_afi(d::AFIInputFile, mesh;
        system = missing,
        phases = missing,
        use_nnc = true,
        active = mesh.cell_map,
        kwarg...
    )
    if ismissing(phases)
        if ismissing(system)
            phases = (AqueousPhase(), LiquidPhase(), VaporPhase())
        else
            phases = JutulDarcy.get_phases(system)
        end
    end
    ncells = number_of_cells(mesh)
    if isnothing(active)
        active = 1:ncells
    end

    data = Dict{String, Vector}()
    IX = d.setup["IX"]
    if haskey(IX, "RESQML")
        resqml = IX["RESQML"]
        if !haskey(resqml, "POROSITY")
            @warn "Porosity is missing in RESQML section of AFI file. Assuming porosity = 1.0 everywhere."
        end
        for (k, resqml_entry) in pairs(resqml)
            if k == "ACTIVE_CELL_FLAG" || k == "GRID"
                continue
            end
            v = vec(resqml_entry["values"])[active]
            set_grid_entry!(data, k, v, ncells)
        end
        found = true
    else
        pillar_grids = find_records(d, "StraightPillarGrid", once = false)
        for pillar_grid in pillar_grids
            found = true
            for (k, v) in pillar_grid.value["CellDoubleProperty"]
                set_grid_entry!(data, k, v, ncells)
            end
        end
    end
    if !found
        @warn "No supported section for grid data found: Did not find RESQML section in AFI file, and no StraightPillarGrid record found. Properties may be missing."
    end
    edits = find_records(d, "BoxPropertyEdit", "IX", steps = false, model = true, once = false)
    if length(edits) > 0
        ijk_lookup = missing
        try
            ijk_lookup = JutulDarcy.ijk_lookup_dict(mesh)
        catch e
            ijk_lookup = Dict{Tuple{Int, Int, Int}, Int}()
            for i in 1:number_of_cells(mesh)
                ijk_lookup[(i, 1, 1)] = i
            end
        end
        for edit in edits
            apply_box_property_edit!(data, edit.value, ncells, ijk_lookup)
        end
    end
    # TODO: Move unit conversion here to properly handle edits that use absolute
    # values and not multipliers.
    domain_kwarg = remap_properties_to_jutuldarcy_names(data, ncells, phases)

    # Heat capacity for components
    d.setup["IX"]
    cmodel = find_records(d, "CompositionalFluidModel", once = true)
    if !isnothing(cmodel)
        lenthalpy = get(cmodel.value, "EnthalpyLiquidHeatCapacity", missing)
        if !ismissing(lenthalpy)
            c1 = get(lenthalpy, "LiquidHeatCapacityCoef1", missing)
            if !ismissing(c1)
                domain_kwarg[:component_heat_capacities] = repeat(c1, 1, ncells)
            end
        end
    end

    if use_nnc
        conn = find_records(d, "ConnectionSet", once = false)
    else
        conn = []
    end
    if length(conn) > 0
        allconn = Dict{String, Any}()
        for c in conn
            nm = c.value["name"]
            if haskey(allconn, nm)
                merge!(allconn[nm], c.value)
            else
                allconn[nm] = c.value
            end
        end
        left = Int[]
        right = Int[]
        trans_nnc = Float64[]
        # global_to_local = afi_to_jutul_cell_map(d, mesh)
        for (k, conn) in pairs(allconn)
            if haskey(conn, "table")
                tab = conn["table"]
                for (c1, c2, t) in zip(tab["Cell1"], tab["Cell2"], tab["Transmissibility"])
                    # TODO: Is this global or local indexing?
                    # push!(left, global_to_local[c1])
                    # push!(right, global_to_local[c2])
                    push!(left, c1)
                    push!(right, c2)
                    push!(trans_nnc, t)
                end
            end
        end
        num_nnc = length(trans_nnc)
        if num_nnc > 0
            @warn "Added $(length(trans_nnc)) NNC connections from AFI file. This functionality is not well tested, especially for meshes with inactive cells."
            domain_kwarg[:nnc] = JutulDarcy.setup_nnc_connections(mesh, left, right, trans_nnc)
        end
    else
        num_nnc = 0
    end

    reservoir = reservoir_domain(mesh; domain_kwarg..., kwarg...)
    # Finally, check if transmissibility override is present
    if haskey(data, "CELL_CENTER_DEPTH")
        cc = reservoir[:cell_centroids, Cells()]
        cc[3, :] .= data["CELL_CENTER_DEPTH"]
    end
    set_transmissibility_override!(reservoir, data, num_nnc)
    return reservoir
end

function set_transmissibility_override!(reservoir, data, num_nnc)
    mesh = physical_representation(reservoir)
    has_trani = haskey(data, "TRANSMISSIBILITY_I")
    has_tranj = haskey(data, "TRANSMISSIBILITY_J")
    has_trank = haskey(data, "TRANSMISSIBILITY_K")
    has_tran_override = has_trani || has_tranj || has_trank
    if has_tran_override
        trans_override = JutulDarcy.reservoir_transmissibility(reservoir)
        ijk = map(cell_ijk(mesh), 1:number_of_cells(mesh))
        if has_trani
            set_trans_override!(trans_override, mesh, ijk, 1, data["TRANSMISSIBILITY_I"], num_nnc)
        end
        if has_tranj
            set_trans_override!(trans_override, mesh, ijk, 2, data["TRANSMISSIBILITY_J"], num_nnc)
        end
        if has_trank
            set_trans_override!(trans_override, mesh, ijk, 3, data["TRANSMISSIBILITY_K"], num_nnc)
        end
        reservoir[:transmissibility_override, Faces()] = trans_override
    end
    return reservoir
end

function set_trans_override!(tran_override, mesh, ijk, dir, tran_xyz, num_nnc)
    N = mesh.faces.neighbors
    active_ix = mesh.cell_map
    if isnothing(active_ix)
        active_ix = 1:number_of_cells(mesh)
    end
    nf = number_of_faces(mesh)
    for (c, val) in enumerate(tran_xyz[active_ix])
        ijk_c = ijk[c][dir]
        if !isfinite(val) || val < 0.0
            continue
        end
        for face in mesh.faces.cells_to_faces[c]
            if face > nf + num_nnc
                continue
            end
            l, r = N[face]
            if l == c
                ijk_other = ijk[r][dir]
            else
                @assert r == c
                ijk_other = ijk[l][dir]
            end
            if ijk_other != ijk_c
                tran_override[face] = val
            end
        end
    end
    return tran_override
end

function apply_box_property_edit!(data, record, ncells, ijk_lookup)
    # @info "??" record
    for (i1, i2, j1, j2, k1, k2, prop, expr) in zip(
            record["I1"],
            record["I2"],
            record["J1"],
            record["J2"],
            record["K1"],
            record["K2"],
            record["Property"],
            record["Expression"],
        )
        (i1 <= i2 && j1 <= j2 && k1 <= k2) || error("Invalid box definition in BoxPropertyEdit record.")
        prop = String(prop)
        v = get_property_from_string(data, prop, ncells; T = Float64)
        apply_box_property_edit_inner!(
            v,
            data,
            i1:i2,
            j1:j2,
            k1:k2,
            prop,
            expr,
            ijk_lookup,
        )
    end
    return data
end

function replace_strings_with_dict_access(s::AbstractString, dictname::AbstractString = "data")
    pat = r"\b([A-Z_]+)\b"
    result = IOBuffer()
    lastidx = 1
    for m in eachmatch(pat, s)
        # Write text before match
        print(result, s[lastidx:m.offset-1])
        # Write replacement
        print(result, "$(dictname)[\"$(m.captures[1])\"][cell_index]")
        lastidx = m.offset + length(m.match)
    end
    # Write any remaining text
    print(result, s[lastidx:end])
    return String(take!(result))
end

function apply_box_property_edit_inner!(vals, data, irange, jrange, krange, prop, expr0, ijk_lookup_dict)
    new_val = tryparse(Float64, expr0)
    if isnothing(new_val)
        # Value is an expression and we have to do eval magic
        expr = replace(expr0, "if" => "ifelse", " " => "")
        expr = replace_strings_with_dict_access(expr, "data")
        F = missing
        try
            F = eval(Meta.parse("($prop, data, cell_index) -> $expr"))
        catch excpt
            @warn "Could not parse expression in BoxPropertyEdit record:\n$expr\nRaw expression: $expr0. Expression will be ignored." excpt
        end
        if !ismissing(F)
            try
                for I in irange, J in jrange, K in krange
                    idx = get(ijk_lookup_dict, (I, J, K), nothing)
                    if !isnothing(idx)
                        new_val = Base.invokelatest(F, vals[idx], data, idx)
                        vals[idx] = new_val
                    end
                end
            catch e
                @warn "Error when applying expression in BoxPropertyEdit record:\n$expr\nRaw expression: $expr0. Expression will be ignored." e
            end
        end
    else
        # The value is a constant number
        for I in irange, J in jrange, K in krange
            idx = get(ijk_lookup_dict, (I, J, K), nothing)
            if !isnothing(idx)
                vals[idx] = new_val
            end
        end
    end
    return vals
end


function set_grid_entry!(data::Dict{String, Vector}, k, v, ncells)
    T = eltype(v)
    if T <: AbstractFloat && T != Float64
        T = Float64
    elseif T <: Integer && T != Int
        T = Int
    end
    nvals = length(v)
    if nvals != ncells
        @warn "Property $k has $nvals values, but the grid has $ncells cells..."
    end
    vals = get_property_from_string(data, k, nvals; T = T)
    vals .= v
    return data
end

function get_property_from_string(data, k, ncells; T = Float64)
    if haskey(data, k)
        v = data[k]
        length(v) == ncells || error("Property $k already exists but has length $(length(v)) instead of $ncells.")
        eltype(v) == T || error("Property $k already exists but has element type $(eltype(v)) instead of $T.")
    else
        if endswith(k, "_MULTIPLIER") 
            v = ones(T, ncells)
        else
            v = zeros(T, ncells)
        end
        data[k] = v
    end
    return v
end

function remap_properties_to_jutuldarcy_names(data, ncells, phases)
    perm = zeros(Float64, 3, ncells)
    poro = ones(Float64, ncells)
    ntg = ones(Float64, ncells)

    function setperm!(ix, vals)
        length(vals) == ncells || error("Permeability component $ix has length $(length(vals)) instead of $ncells.")
        perm[ix, :] .= vals
    end
    out = Dict{Symbol, Any}(
        :permeability => perm,
        :porosity => poro,
        :net_to_gross => ntg,
    )
    for (k, vals) in pairs(data)
        # TODO: pore volumes, transmissibilities, etc...
        if k == "PERM_I" || k == "PERM_X"
            setperm!(1, vals)
        elseif k == "PERM_J" || k == "PERM_Y"
            setperm!(2, vals)
        elseif k == "PERM_K" || k == "PERM_Z"
            setperm!(3, vals)
        elseif k == "POROSITY"
            poro .= vals
        elseif k == "NET_TO_GROSS_RATIO"
            ntg .= vals
        elseif k == "PORE_VOLUME"
            out[:pore_volume] = vals
        elseif k == "PORE_VOLUME_MULTIPLIER"
            out[:pore_volume_multiplier] = vals
        elseif k in ["TRANSMISSIBILITY_I", "TRANSMISSIBILITY_J", "TRANSMISSIBILITY_K"]
            # Handled separately
            continue
        elseif k == "THERMAL_CONDUCTIVITY_ROCK"
            out[:rock_thermal_conductivity] = vals
        elseif k == "ROCK_HEAT_CAPACITY"
            out[:rock_heat_capacity] = vals
        elseif startswith(k, "THERMAL_CONDUCTIVITY_")
            if !haskey(out, :fluid_thermal_conductivity)
                out[:fluid_thermal_conductivity] = handle_suffixed_entries(data, ncells, "THERMAL_CONDUCTIVITY_", phases)
            end
        else
            out[Symbol(k)] = vals
        end
    end
    return out
end

function handle_suffixed_entries(data, nc, prefix, phases = (AqueousPhase(), LiquidPhase(), VaporPhase()))
    nph = length(phases)
    val = zeros(nph, nc)
    for (phno, phase) in enumerate(phases)
        if phase == AqueousPhase()
            suffix = "WATER"
        elseif phase == LiquidPhase()
            suffix = "OIL"
        elseif phase == VaporPhase()
            suffix = "GAS"
        else
            error("Unsupported phase type $phase for suffixed property handling.")
        end
        val[phno, :] = data["$prefix$suffix"]
    end
    return val
end
