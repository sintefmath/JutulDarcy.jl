# Keep geometry by reference while isolating DFM's volume and property updates.
function copy_fracture_data(domain::DataDomain)
    return DataDomain(physical_representation(domain), copy(domain.entities), deepcopy(domain.data))
end

function setup_unified_dfm_model(matrix, fractures, system;
        wells = [], extend_wells = true, subtract_fracture_volume = true,
        extra_out = false, kwarg...)
    # Initialize the same coefficients and well completions as coupled DFM.
    coupled = setup_fractured_reservoir_model(copy_fracture_data(matrix),
        copy_fracture_data(fractures), system;
        wells = [copy_fracture_data(well) for well in wells],
        extend_wells = extend_wells, subtract_fracture_volume = subtract_fracture_volume,
        kwarg...)
    parameters = setup_parameters(coupled)
    domain = unified_dfm_domain(coupled, parameters)
    nm = number_of_cells(matrix)
    wells = [unified_fracture_well(well, parameters[name], nm)
        for (name, well) in pairs(get_model_wells(coupled; data_domain = true))]
    model = setup_reservoir_model(domain, system; wells = wells, kwarg...)
    if JutulDarcy.model_is_thermal(reservoir_model(model))
        set_parameters!(reservoir_model(model), FluidThermalConductivities = EFMFluidThermalConductivities())
    end
    out = model
    if extra_out
        out = (model, setup_parameters(model))
    end
    return out
end

function unified_dfm_domain(model, parameters)
    matrix = reservoir_domain(model)
    fractures = model[:Fractures].data_domain
    matrix_mesh = physical_representation(matrix)
    fracture_mesh = physical_representation(fractures)
    nm = number_of_cells(matrix)
    fmc = FractureMatrixConnection()
    mc, fc = fractures[:matrix_cells, fmc], fractures[:connection_cells, fmc]
    mf = fractures[:matrix_faces, fmc]
    N = hcat(matrix[:neighbors], fractures[:neighbors] .+ nm, hcat(mc, fc .+ nm)')
    volumes = vcat(matrix[:volumes], fractures[:volumes])
    centers = hcat(matrix[:cell_centroids], fractures[:cell_centroids])
    connection_centers = (matrix[:cell_centroids][:, mc]+fractures[:cell_centroids][:, fc])/2
    face_centers = hcat(matrix[:face_centroids], fractures[:face_centroids], connection_centers)
    normals = hcat(matrix[:normals], fractures[:normals], matrix[:normals][:, mf])
    areas = vcat(matrix[:areas], fractures[:areas], tpfv_geometry(fracture_mesh).volumes[fc])
    mg = tpfv_geometry(matrix_mesh)
    geo = TwoPointFiniteVolumeGeometry(N, areas, volumes, normals, centers, face_centers;
        boundary_areas = mg.boundary_areas, boundary_normals = mg.boundary_normals,
        boundary_centroids = mg.boundary_centroids, boundary_neighbors = mg.boundary_neighbors)
    mesh = EmbeddedReservoirMesh(matrix_mesh, fracture_mesh, geo, geometry_tags(geo))
    domain = reservoir_domain(mesh)
    parents = [fracture_parent_cell(fracture_mesh, matrix_mesh, cell)
        for cell in 1:number_of_cells(fractures)]
    for (key, (value, entity)) in pairs(matrix)
        if entity == Cells()
            if haskey(fractures, key, Cells())
                fracture_values = fractures[key, Cells()]
            elseif key == :net_to_gross
                fracture_values = ones(eltype(value), length(parents))
            else
                fracture_values = fracture_cell_subset(value, parents)
            end
            domain[key, entity] = append_cell_property(value, fracture_values)
        elseif entity == BoundaryFaces() || entity == Jutul.NoEntity()
            domain[key, entity] = deepcopy(value)
        end
    end
    domain[:volumes, Cells()] = volumes
    domain[:cell_centroids, Cells()] = centers
    domain[:is_fracture, Cells()] = vcat(zeros(Int, nm), ones(Int, number_of_cells(fractures)))
    if haskey(fractures, :fracture_id, Cells())
        domain[:fracture_id, Cells()] = vcat(zeros(Int, nm), fractures[:fracture_id])
    end
    pm, pf = parameters[:Reservoir], parameters[:Fractures]
    domain[:transmissibilities, Faces()] = vcat(pm[:Transmissibilities],
        pf[:Transmissibilities], pf[:FractureMatrixTransmissibility])
    if JutulDarcy.model_is_thermal(reservoir_model(model))
        # DFM cross terms use total heat conductivity independent of saturation.
        # Store it in the rock term, with zero fluid conductivity on these faces.
        domain[:rock_thermal_conductivities, Faces()] = vcat(pm[:RockThermalConductivities],
            pf[:RockThermalConductivities], pf[:FractureMatrixThermalConductivity])
        nph = size(pm[:FluidThermalConductivities], 1)
        domain[:fluid_thermal_conductivities, Faces()] = hcat(pm[:FluidThermalConductivities],
            pf[:FluidThermalConductivities], zeros(eltype(volumes), nph, length(mc)))
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
