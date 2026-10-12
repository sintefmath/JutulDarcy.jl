# Geometry/property assembly shared by EFM and DFM. Their flow and heat
# coefficients are initialized separately to preserve each discretization.
function unified_fracture_domain(matrix, fractures, connection_areas, connection_normals;
        volumes = vcat(matrix[:volumes], fractures[:volumes]))
    matrix_mesh = physical_representation(matrix)
    fracture_mesh = physical_representation(fractures)
    nm, nf = number_of_cells(matrix), number_of_cells(fractures)
    fmc = FractureMatrixConnection()
    mc, fc = fractures[:matrix_cells, fmc], fractures[:connection_cells, fmc]
    N = hcat(matrix[:neighbors], fractures[:neighbors] .+ nm, hcat(mc, fc .+ nm)')
    centers = hcat(matrix[:cell_centroids], fractures[:cell_centroids])
    connection_centers = (matrix[:cell_centroids][:, mc]+fractures[:cell_centroids][:, fc])/2
    face_centers = hcat(matrix[:face_centroids], fractures[:face_centroids], connection_centers)
    normals = hcat(matrix[:normals], fractures[:normals], connection_normals)
    areas = vcat(matrix[:areas], fractures[:areas], connection_areas)
    mg = tpfv_geometry(matrix_mesh)
    geo = TwoPointFiniteVolumeGeometry(N, areas, volumes, normals, centers, face_centers;
        boundary_areas = mg.boundary_areas, boundary_normals = mg.boundary_normals,
        boundary_centroids = mg.boundary_centroids, boundary_neighbors = mg.boundary_neighbors)
    mesh = EmbeddedReservoirMesh(matrix_mesh, fracture_mesh, geo, geometry_tags(geo))
    domain = reservoir_domain(mesh)
    parents = [fracture_parent_cell(fracture_mesh, matrix_mesh, cell) for cell in 1:nf]
    for (key, (value, entity)) in pairs(matrix)
        if entity == Cells()
            if haskey(fractures, key, Cells())
                fracture_values = fractures[key, Cells()]
            elseif key == :net_to_gross
                fracture_values = ones(eltype(value), nf)
            elseif value isa AbstractVector || value isa AbstractMatrix
                fracture_values = fracture_cell_subset(value, parents)
            else
                throw(ArgumentError("cannot extend matrix cell property $key with type $(typeof(value)); supply the corresponding fracture property"))
            end
            domain[key, entity] = append_cell_property(value, fracture_values)
        elseif entity == BoundaryFaces() || entity == Jutul.NoEntity()
            domain[key, entity] = deepcopy(value)
        end
    end
    domain[:volumes, Cells()] = volumes
    domain[:cell_centroids, Cells()] = centers
    domain[:is_fracture, Cells()] = vcat(zeros(Int, nm), ones(Int, nf))
    if haskey(fractures, :fracture_id, Cells())
        domain[:fracture_id, Cells()] = vcat(zeros(Int, nm), fractures[:fracture_id])
    end
    return domain
end

function unified_fracture_well(well::DataDomain, parameters, nm)
    original = physical_representation(well)
    unified = well
    if haskey(original.perforations, :fracture)
        w = deepcopy(original)
        perforations = (self = vcat(w.perforations.self, w.perforations.self_fracture),
            reservoir = vcat(w.perforations.reservoir, w.perforations.fracture .+ nm))
        representation = JutulDarcy.MultiSegmentWell(w.type, w.num_nodes, w.num_segments,
            length(perforations.self), perforations, w.neighborship, w.end_nodes,
            w.surface, w.name, w.segment_models, w.include_inertia, w.multiwell)
        unified = DataDomain(representation)
        p, fp = JutulDarcy.Perforations(), JutulDarcy.FracturePerforations()
        for (key, (value, entity)) in pairs(well)
            if entity == p
                fracture_values = well[Symbol(key, :_frac), fp]
                unified[key, p] = append_cell_property(value, fracture_values)
            elseif entity != fp
                unified[key, entity] = deepcopy(value)
            end
        end
        unified[:well_index, p] = vcat(parameters[:WellIndices], parameters[:FractureWellIndices])
        if haskey(parameters, :WellIndicesThermal)
            unified[:thermal_well_index, p] = vcat(parameters[:WellIndicesThermal], parameters[:FractureWellIndicesThermal])
        end
    end
    return unified
end

function unified_fracture_well(well::DataDomain, nm)
    parameters = Dict{Symbol, Any}()
    if haskey(physical_representation(well).perforations, :fracture)
        parameters[:WellIndices] = Jutul.default_parameter_values(well, nothing,
            JutulDarcy.WellIndices(), :WellIndices)
        parameters[:FractureWellIndices] = Jutul.default_parameter_values(well, nothing,
            JutulDarcy.FractureWellIndices(), :FractureWellIndices)
        parameters[:WellIndicesThermal] = Jutul.default_parameter_values(well, nothing,
            JutulDarcy.WellIndicesThermal(), :WellIndicesThermal)
        parameters[:FractureWellIndicesThermal] = Jutul.default_parameter_values(well, nothing,
            JutulDarcy.FractureWellIndicesThermal(), :FractureWellIndicesThermal)
    end
    unified = unified_fracture_well(well, parameters, nm)
    return unified
end

function fracture_parent_cell(fractures::EmbeddedFractureMesh, matrix, cell)
    return first(fractures.geometry.patches[cell].matrix_cells)
end

function fracture_parent_cell(fractures::Jutul.EmbeddedMesh, matrix, cell)
    if cell > length(fractures.parent_faces)
        intersection = findfirst(==(cell), fractures.intersection_cells)
        cell = first(fractures.intersection_neighbors[intersection])
    end
    face = fractures.parent_faces[cell]
    return get_neighborship(matrix)[1, face]
end
