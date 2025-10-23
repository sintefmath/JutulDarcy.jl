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
        use_nnc = true,
        active = mesh.cell_map
    )
    ncells = number_of_cells(mesh)
    if isnothing(active)
        active = 1:ncells
    end
    # perm = zeros(Float64, 3, ncells)
    # poro = ones(Float64, ncells)
    # ntg = ones(Float64, ncells)

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
    # TODO: Move unit conversion here to properly handle edits.
    if use_nnc
        conn = find_records(d, "ConnectionSet", once = false)
    else
        conn = []
    end
    has_trani = haskey(data, "TRANSMISSIBILITY_I")
    has_tranj = haskey(data, "TRANSMISSIBILITY_J")
    has_trank = haskey(data, "TRANSMISSIBILITY_K")
    if length(conn) > 0 || has_trani || has_tranj || has_trank
        # NNC
        # TRANSMIBBILITY_(I/J/K)
        error("Not finished.")
    end
    domain_kwarg = remap_properties_to_jutuldarcy_names(data, ncells)
    return reservoir_domain(mesh; domain_kwarg...)
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

function remap_properties_to_jutuldarcy_names(data, ncells)
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
        else
            out[Symbol(k)] = vals
        end
    end
    return out
end
