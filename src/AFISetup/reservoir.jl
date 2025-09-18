function JutulDarcy.reservoir_domain(d::AFIInputFile; mesh = missing)
    mesh = setup_mesh_afi(d, mesh)

    # TODO: Handle regions and better warnings
    return setup_reservoir_domain_afi(d, mesh)
end

function set_grid_entry!(data, k, v)
    perm = data[:permeability]
    poro = data[:porosity]
    ntg = data[:net_to_gross]
    if k == "PERM_I" || k == "PERM_X"
        perm[1, :] .= v
    elseif k == "PERM_J" || k == "PERM_Y"
        perm[2, :] .= v
    elseif k == "PERM_K" || k == "PERM_Z"
        perm[3, :] .= v
    elseif k == "POROSITY"
        poro .= v
    elseif k == "NET_TO_GROSS_RATIO"
        ntg .= v
    else
        data[Symbol(k)] = v
    end
    return data
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

function setup_mesh_afi(afi::AFIInputFile, mesh)
    IX = afi.setup["IX"]
    if haskey(IX, "RESQML")
        resqml = IX["RESQML"]
        if ismissing(mesh)
            haskey(IX["RESQML"], "GRID") || error("No GRID section in converted AFI file under the RESQML section.")
            mesh = GeoEnergyIO.mesh_from_grid_section(resqml["GRID"])
        end
    else
        pillar_grid = find_records(afi, "StraightPillarGrid", once = true)
        if isnothing(pillar_grid)
            error("No RESQML section in AFI file, and no StraightPillarGrid record found.")
        end
        cartdims = cartdims_from_structured_info(afi)
        nx, ny, nz = cartdims
        dx = pillar_grid.value["DeltaX"]
        dy = pillar_grid.value["DeltaY"]
        dz = pillar_grid.value["DeltaZ"]
        tops = pillar_grid.value["PillarTops"]
        grid = Dict(
            "cartDims" => cartdims,
            "DXV" => dx,
            "DYV" => dy,
            "DZV" => dz,
            "TOPS" => tops[1:nx*ny],
        )
        mesh = mesh_from_grid_section(grid)
    end
    !ismissing(mesh) || error("Could not set up mesh from AFI file.")
    return mesh
end

function setup_reservoir_domain_afi(d::AFIInputFile, mesh)
    active = mesh.cell_map
    ncells = number_of_cells(mesh)
    perm = zeros(Float64, 3, ncells)
    poro = ones(Float64, ncells)
    ntg = ones(Float64, ncells)
    domain_kwarg = Dict{Symbol, Any}(
        :permeability => perm,
        :porosity => poro,
        :net_to_gross => ntg,
    )

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
            set_grid_entry!(domain_kwarg, k, v)
        end
        found = true
    else
        pillar_grid = find_records(afi, "StraightPillarGrid", once = true)
        if !isnothing(pillar_grid)
            found = true
            for (k, v) in pillar_grid.value["CellDoubleProperty"]
                set_grid_entry!(domain_kwarg, k, v)
            end
        end
    end
    if !found
        @warn "No supported section for grid data found: Did not find RESQML section in AFI file, and no StraightPillarGrid record found. Properties may be missing."
    end
    return reservoir_domain(mesh; domain_kwarg...)
end
