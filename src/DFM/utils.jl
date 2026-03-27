# using LinearAlgebra
# using StaticArrays

function fracture_domain(mesh::JutulMesh, matrix::DataDomain;
    connection_cells_fracture=missing,
    matrix_faces=missing,
    matrix_cells=missing,
    kwarg...)

    # Define fracture cells connected to matrix (default: all except intersection cells)
    if ismissing(connection_cells_fracture)
        connection_cells_fracture = 1:number_of_cells(mesh)
        if hasproperty(mesh, :intersection_cells)
            connection_cells_fracture = setdiff(
                connection_cells_fracture, mesh.intersection_cells)
        end
    end
    # Define matrix faces corresponding to fractures (if given)
    if ismissing(matrix_faces)
        if hasproperty(mesh, :parent_faces)
            matrix_faces = mesh.parent_faces
        else
            @warn "Matrix faces not provided and not found in fracture mesh. \
            Fractures will not affect flow through matrix cells unless already \
            accounted for in the provided matrix domain."
        end
    end
    # Define matrix cells connected to fractures (if given)
    if ismissing(matrix_cells)
        !ismissing(matrix_faces) || error("Must provide matrix_cells if \
            matrix_faces is not provided in order to set up matrix/fracture \
            connections.")
        matrix_mesh = physical_representation(matrix)
        N = get_neighborship(matrix_mesh)
        matrix_cells = N[:, matrix_faces][:]
        matrix_faces = repeat(matrix_faces, inner=2)
        connection_cells_fracture = repeat(connection_cells_fracture, inner=2)
    end
    length(connection_cells_fracture) == length(matrix_faces) == length(matrix_cells) ||
    error("Length of connection_cells_fracture, matrix_faces, and matrix_cells must match")

    # Set up fracture domain
    fracture = fracture_domain(mesh; kwarg...)
    # Add fracture-matrix connection entitiy to fracture domain
    fmc = FractureMatrixConnection()
    fmc = JutulDarcy.FractureMatrixConnection()
    num_fmc = length(matrix_faces)
    fracture.entities[fmc] = num_fmc

    # Add matrix properties to fracture domain for use in cross-term calculations
    x = matrix[:cell_centroids]
    K = matrix[:permeability]
    ϕ = matrix[:porosity]
    Λ_f = matrix[:fluid_thermal_conductivity]
    Λ_r = matrix[:rock_thermal_conductivity]
    Λ_r = vec(Λ_r)
    ϕ = vec(ϕ)
    fracture[:connection_cells, fmc] = connection_cells_fracture
    fracture[:matrix_faces, fmc] = matrix_faces
    fracture[:matrix_cells, fmc] = matrix_cells
    fracture[:matrix_cell_centroids, fmc] = x[:, matrix_cells]
    fracture[:matrix_permeability, fmc] = K[matrix_cells]
    fracture[:matrix_fluid_thermal_conductivity, fmc] = Λ_f[matrix_cells]
    fracture[:matrix_rock_thermal_conductivity, fmc] = Λ_r[matrix_cells]
    fracture[:matrix_porosity, fmc] = ϕ[matrix_cells]
    cell_normals = fill(NaN, size(x, 1), number_of_cells(mesh))
    for cell in connection_cells_fracture
        cell_normals[:, cell] = Jutul.EmbeddedMeshes.cell_normal(mesh, cell)
    end
    fracture[:cell_normals, Cells()] = cell_normals

    return fracture

end

function fracture_domain(mesh::JutulMesh;
    aperture=0.5e-3si_unit(:meter),
    hydralic_aperture=aperture,
    permeability=missing,
    kwarg...
    )

    all(isfinite, aperture) || throw(ArgumentError(
        "Keyword argument aperture has non-finite entries."))
    minimum(aperture) >= 0 || throw(ArgumentError(
        "All aperture values must be non-negative."))
    all(isfinite, hydralic_aperture) || throw(ArgumentError(
        "Keyword argument hydralic_aperture has non-finite entries."))
    minimum(hydralic_aperture) >= 0 || throw(ArgumentError(
        "All hydralic_aperture values must be non-negative."))

    if ismissing(permeability)
        # Use cubic law to compute permeability from hydraulic aperture if not provided
        permeability = (hydralic_aperture.^2)./12.0
    end
    # Set up the fracture domain with the provided aperture and computed permeability
    domain = reservoir_domain(mesh;
        permeability = permeability,
        kwarg...
    )
    # Set aperture
    domain[:aperture, Cells()] = aperture
    # Adjust cell volumes and face areas to account for aperture
    geo = tpfv_geometry(mesh)
    aperture = domain[:aperture, Cells()]
    volumes = geo.volumes.*aperture
    volumes[mesh.intersection_cells] .*= aperture[mesh.intersection_cells]
    domain[:volumes, Cells()] = volumes
    N = domain[:neighbors]
    face_aperture = vec(sum(aperture[N], dims=1)./2)
    domain[:areas, Faces()] = geo.areas.*face_aperture

    return domain

end

function setup_fractured_reservoir_model(matrix::DataDomain, fractures::DataDomain, system::Union{JutulSystem, Symbol};
    wells=[], kwarg...)

    fmc = JutulDarcy.FractureMatrixConnection()
    if haskey(fractures, :matrix_faces)
        matrix_faces = fractures[:matrix_faces, fmc]
    else
        @warn "Matrix faces corresponding to fractures not found in fracture \
        domain. Fractures will not affect flow accross any matrix faces unless \
        already accounted for in the matrix domain."
        matrix_faces = missing
    end
    if haskey(fractures, :matrix_cells)
        matrix_cells = fractures[:matrix_cells, fmc]
    else
        @warn "Matrix cells connected to fractures not found in fracture \
        domain. Cannot set up matrix/fracture cross terms and adjust matrix \
        cell volumes to account for fracture volumes."
        matrix_cells = missing
    end

    if ismissing(matrix_faces)        
        fractured_wells = wells
    else
        matrix_mesh = physical_representation(matrix)
        fractured_wells = []
        for well in wells
            well = add_fractures_to_well(well, fractures, matrix_mesh, unique(matrix_faces))
            push!(fractured_wells, well)
        end
    end

    # Get matrix model
    model = setup_reservoir_model(matrix, system; wells=fractured_wells, kwarg...)
    matrix_model = model.models[:Reservoir]

    # Set up fracture model with same system and context as matrix model
    fracture_model = setup_reservoir_model(fractures, system;
        context=model.context, kwarg...)
    if fracture_model isa Jutul.MultiModel
        fracture_model = fracture_model.models[:Reservoir]
    end
    set_parameters!(fracture_model, Transmissibilities = TransmissibilitiesDFM())
    thermal = model_is_thermal(matrix_model)
    if thermal
        set_parameters!(fracture_model,
            RockThermalConductivities = RockThermalConductivitiesDFM(),
            FluidThermalConductivities = FluidThermalConductivitiesDFM(),
        )
    end
    # Add extra entities to fracture model for fracture-matrix connections
    fracture_model.domain.entities[fmc] = count_entities(fractures, fmc)

    # Block flow across matrix faces that are connected to fractures
    if !ismissing(matrix_faces)
        block_fracture_face_connections!(matrix_model, unique(matrix_faces))
    end
    if !ismissing(matrix_cells)
        adjust_matrix_cell_volumes!(matrix_model, fractures)
    end

    # Add the fracture model to the multimodel
    model = add_fracture_model(model, fracture_model)

    # Set up DFM cross-terms
    add_fracture_cross_terms!(model, fractures)

    return model

end

function block_fracture_face_connections!(matrix_model::SimulationModel, faces)

    thermal = model_is_thermal(matrix_model)
    matrix = matrix_model.data_domain

    # Set zero transmissibility for fracture faces
    if haskey(matrix, :transmissibilities)
        T = matrix[:transmissibilities, Faces()]
    else
        T = Jutul.default_parameter_values(
            matrix, matrix_model, Transmissibilities(), :Transmissibilities)
    end
    T[faces] .= 0.0
    matrix[:transmissibilities, Faces()] = T
    
    thermal || return

    # Set zero rock thermal conductivity for fracture faces
    if haskey(matrix, :rock_thermal_conductivities)
        Λr = matrix[:rock_thermal_conductivities, Faces()]
    else
        Λr = Jutul.default_parameter_values(
            matrix, matrix_model, RockThermalConductivities(), :RockThermalConductivities)
    end
    Λr[faces] .= 0.0
    matrix[:rock_thermal_conductivities, Faces()] = Λr

    # Set zero fluid thermal conductivity for fracture faces
    if haskey(matrix, :fluid_thermal_conductivities)
        Λf = matrix[:fluid_thermal_conductivities, Faces()]
    else
        Λf = Jutul.default_parameter_values(
            matrix, matrix_model, FluidThermalConductivities(), :FluidThermalConductivities)
    end
    Λf[faces] .= 0.0
    matrix[:fluid_thermal_conductivities, Faces()] = Λf

    return

end

function adjust_matrix_cell_volumes!(matrix_model::SimulationModel, fractures::DataDomain)

    matrix = matrix_model.data_domain

    fracture_volumes = fractures[:volumes, Cells()]
    matrix_volumes = matrix[:volumes, Cells()]
    fmc = JutulDarcy.FractureMatrixConnection()
    fracture_cells = fractures[:connection_cells, fmc]
    matrix_cells = fractures[:matrix_cells, fmc]
    for (cf, mc) in zip(fracture_cells, matrix_cells)
        vf = fracture_volumes[cf]
        matrix_volumes[mc] -= vf/2
    end
    matrix[:volumes, Cells()] = matrix_volumes

    return

end

function add_fractures_to_well(well, fractures::DataDomain, matrix_mesh, matrix_faces; kwarg...)

    well0 = well
    neighborship, perforation_cells_fractures, perforation_cells_fractures_self, old_ix = 
    add_fracture_cells_to_well(matrix_mesh, well0.representation, matrix_faces)

    has_fracture_perforations = length(perforation_cells_fractures) > 0
    if !has_fracture_perforations
        return well0
    end

    num_segments = size(neighborship, 2)
    segment_models = [SegmentWellBoreFrictionHB() for _ in 1:num_segments]
    for seg = 1:num_segments
        if old_ix[seg] != 0
            segment_models[seg] = well0.representation.segment_models[old_ix[seg]]
        end
    end

    perf = (
        self = well0.representation.perforations.self,
        reservoir = well0.representation.perforations.reservoir,
        self_fracture = perforation_cells_fractures_self,
        fracture = perforation_cells_fractures,
    )

    # TODO: Find radii and direction nearest fracture perforations
    radius = well0[:perforation_radius, Perforations()][1]
    dir = well0[:perforation_direction, Perforations()][1]
    K_fractures = fractures[:permeability]
    # Compute effective thermal conductivity
    Λf = fractures[:fluid_thermal_conductivity]
    Λr = fractures[:rock_thermal_conductivity]
    if haskey(fractures, :net_to_gross)
        ntg = fractures[:net_to_gross]
    else
        ntg = missing
    end
    ϕ = fractures[:porosity]
    Λr = vec(Λr)
    ϕ = vec(ϕ)
    if size(Λf, 1) == 1 || Λf isa Vector
        Λf = vec(Λf)
        Λ = ϕ.*Λf + (1.0 .- ϕ).*Λr
    else
        # TODO: This is a bit of a hack. We should really have a proper way to
        # do this inside the equations for multiphase flow.
        Λ = Λr
    end

    if has_fracture_perforations
        fracture_mesh = physical_representation(fractures)
        tmp_frac = setup_well(fracture_mesh, K_fractures, perforation_cells_fractures;
        simple_well=false, radius=radius, dir=dir, net_to_gross=ntg,
        thermal_conductivity=Λ,
        kwarg...)
    else
        tmp_frac = missing
    end

    W = MultiSegmentWell(
        well0.representation.type,
        maximum(neighborship, init = 1),
        length(segment_models),
        length(well0.representation.perforations.self),
        perf,
        neighborship,
        well0.representation.end_nodes,
        well0.representation.surface,
        well0.representation.name,
        segment_models,
    )

    well = DataDomain(W)
    for (k, v) in well0.data
        val, entity = v
        if ismissing(tmp_frac)
            well[k, entity] = val
            continue
        end
        frac_val = tmp_frac[k, entity]
        if entity == Cells()
            # Cell data
            if k ∈ [:cell_centroids, :perforation_centroids]
                new_val = hcat(val, frac_val)
            else
                new_val = vcat(val, frac_val)
            end
        elseif entity == Faces()
            new_val = similar(val, count_entities(W, Faces()))
            new_val[old_ix .!= 0] = val[old_ix[old_ix .!= 0]]
            fixed = old_ix .!= 0
            # Insert nearby values for new fracture perforation faces
            while any(.!fixed)
                for seg in 1:num_segments
                    if !fixed[seg]
                        pseg = max(seg-1, 1)
                        nseg = min(seg+1, num_segments)
                        if fixed[pseg]
                            new_val[seg] = new_val[pseg]
                            fixed[seg] = true
                        elseif fixed[nseg]
                            new_val[seg] = new_val[nseg]
                            fixed[seg] = true
                        end
                    end
                end
            end
        elseif entity == Perforations()
            well[k, entity] = val
            new_val = frac_val
            k = Symbol(String(k)*"_frac")
            entity = FracturePerforations()
        end
        well[k, entity] = new_val
    end

    wc = well.representation.perforations.self_fracture
    fc = well.representation.perforations.fracture
    
    well[:cell_length][wc] = fractures[:aperture][fc]
    Δ = well[:cell_dims_frac, FracturePerforations()]
    dir = well[:perforation_direction_frac, FracturePerforations()]
    for (c, (Δ_k, dir_k)) in enumerate(zip(Δ, dir))
        if dir_k isa Symbol
            d = findfirst(==(dir_k), [:x, :y, :z])
        else
            d = last(findmax(abs.(dir_k)))
        end
        Δ_k = [Δ_k...]
        Δ_k[d] = fractures[:aperture][fc[c]]
        Δ[c] = Tuple(Δ_k)
    end
    well[:cell_dims_frac, FracturePerforations()] = Δ

    return well

end

function add_fracture_cells_to_well(g_matrix, representation, fracture_faces)
    
    N = representation.neighborship
    nc = maximum(N)
    w2m = fill(0, nc)
    perf = representation.perforations
    w2m[perf.self] .= perf.reservoir

    N_new = Matrix{Int}(undef, 2, 0)
    N_matrix = get_neighborship(g_matrix)
    perforation_cells_fractures = Int[]
    perforation_cells_fractures_self = Int[]

    wc = maximum(N)
    added = false
    old_ix = Int[]
    for (seg, mcells) in enumerate(eachcol(w2m[N]))
        ix = vec(all(mcells .== N_matrix, dims = 1) .|| all(reverse(mcells) .== N_matrix, dims = 1))
        if any(ix)
            @assert sum(ix) == 1 "Failed to find unique matching face for well segment perforation"
            face = findfirst(ix)
            if face ∈ fracture_faces
                wc += 1
                fcell = findfirst(fracture_faces .== face)
                push!(perforation_cells_fractures, fcell)
                push!(perforation_cells_fractures_self, wc)
                m2f = [N[1, seg], wc]
                f2m = [wc, N[2, seg]]
                N_new = hcat(N_new, m2f, f2m)
                push!(old_ix, 0, 0)
                added = true
            end
        end
        if !added
            N_new = hcat(N_new, N[:, seg])
            push!(old_ix, seg)
        end
        added = false
    end

    return N_new, perforation_cells_fractures, perforation_cells_fractures_self, old_ix

end

function add_fracture_model(model::MultiModel, fracture_model::SimulationModel)

    new_models = JutulStorage()
    new_models[:Reservoir] = model.models[:Reservoir]
    new_models[:Fractures] = fracture_model
    groups = Int[1, 1]
    for (k, v) in pairs(model.models)
        if k != :Reservoir
            new_models[k] = v
            push!(groups, 2)
        end
    end
    if isnothing(model.groups)
        groups = nothing
    end
    model0 = model
    model = Jutul.MultiModel(new_models; groups=groups)
    for ct in model0.cross_terms
        push!(model.cross_terms, ct)
    end

    return model

end

function add_fracture_cross_terms!(model::MultiModel, fractures::DataDomain)
    # Set up DFM cross-terms

    fmc = JutulDarcy.FractureMatrixConnection()

    set_parameters!(model.models[:Fractures],
        FractureMatrixTransmissibility = FractureMatrixTransmissibility(),
        FractureMatrixGravityDifference = FractureMatrixGravityDifference(),
    )
    matrix_cells = fractures[:matrix_cells, fmc]
    fracture_cells = fractures[:connection_cells, fmc]
    ct = MatrixFromFractureFlowCT(matrix_cells, fracture_cells)
    add_cross_term!(model, ct, target = :Reservoir, source = :Fractures, equation = :mass_conservation)

    thermal = model_is_thermal(model.models[:Reservoir])
    if thermal
        set_parameters!(model.models[:Fractures],
            FractureMatrixThermalConductivity = FractureMatrixThermalConductivity(),
        )
        ct = MatrixFromFractureThermalCT(matrix_cells, fracture_cells)
        add_cross_term!(model, ct, target = :Reservoir, source = :Fractures, equation = :energy_conservation)
    end

    for (name, well_model) in get_model_wells(model)
        g = physical_representation(well_model)
        if g isa WellDomain
            if !haskey(g.perforations, :fracture)
                continue
            end
            fc = vec(g.perforations.fracture)
            wc = vec(g.perforations.self_fracture)
            set_parameters!(model.models[name],
                FractureWellIndices = FractureWellIndices()
            )
            ct = FracturesFromWellFlowCT(fc, wc)
            add_cross_term!(model, ct, target = :Fractures, source = name, equation = :mass_conservation)
            if thermal
                set_parameters!(model.models[name],
                    FractureWellIndicesThermal = FractureWellIndicesThermal()
                )
                ct = FracturesFromWellThermalCT(fc, wc)
                add_cross_term!(model, ct, target = :Fractures, source = name, equation = :energy_conservation)
            end
        end
    end

end

function reservoir_conductivity_dfm(domain, Λ)

    mesh = physical_representation(domain)
    star_delta = isempty(mesh.intersection_cells) # Use star-delta if intersection cells are eliminated
    N = domain[:neighbors]
    nc = number_of_cells(mesh)
    faces, facepos = get_facepos(N, nc)
    T_hf = compute_half_face_trans(mesh, 
        domain[:cell_centroids, Cells()],
        domain[:face_centroids, Faces()], 
        domain[:areas, Faces()],
        Λ,
        domain[:aperture, Cells()],
        faces, facepos)
    T = Jutul.EmbeddedMeshes.compute_face_trans_dfm(
            T_hf, N, mesh.intersection_neighbors, star_delta)

    return T

end

function setup_fractured_reservoir_model___(matrix::DataDomain, fractures::DataDomain, system::Union{JutulSystem, Symbol};
    block_backend = true,
    kwarg...)

    matrix_faces = unique(fractures[:matrix_faces])
    # Set zero transmissibility accross matrix cells that are fractures
    T = reservoir_transmissibility(matrix)
    T[matrix_faces] .= 0.0
    matrix[:transmissibilities, Faces()] = T

    if haskey(matrix, :rock_thermal_conductivity)
        T = reservoir_conductivity(matrix)
        T[matrix_faces] .= 0.0
        matrix[:rock_thermal_conductivities, Faces()] = T
    end

    if haskey(matrix, :fluid_thermal_conductivity)
        nph = 1
        @warn "Assuming single-phase system for setting up fluid thermal conductivity"
        C = matrix[:fluid_thermal_conductivity]
        phi = matrix[:porosity]
        if C isa Vector
            T = compute_face_trans(matrix, phi.*C)
            T = repeat(T', nph, 1)
        else
            size(C, 1) == nph || error("Expected size $(nph) x num_cells for :fluid_thermal_conductivity, got size $(size(C))")
            nf = number_of_faces(matrix)
            T = zeros(nph, nf)
            for ph in 1:nph
                T[ph, :] = compute_face_trans(matrix, phi.*C[ph, :])
            end
        end
        T[matrix_faces] .= 0.0
        matrix[:fluid_thermal_conductivities, Faces()] = T
    end

    # Adjust matrix cell volumes to account for fracture volumes
    fracture_volumes = fractures[:volumes, Cells()]
    matrix_volumes = matrix[:volumes, Cells()]
    N = get_neighborship(physical_representation(matrix))
    for (cf, face) in enumerate(matrix_faces)
        vf = fracture_volumes[cf]
        matrix_volumes[N[:, face]] .-= vf/2
    end
    matrix[:volumes, Cells()] = matrix_volumes
    
    model = setup_reservoir_model(matrix, system; block_backend = block_backend, kwarg...)
    has_thermal = haskey(model[:Reservoir].equations, :energy_conservation)

    if has_thermal
        fmesh = physical_representation(fractures)
        geo = tpfv_geometry(fmesh)
        nc = number_of_cells(fmesh)
        N = get_neighborship(fmesh)
        faces, facepos = get_facepos(N, nc)
        ϕ = fractures[:porosity]
        for (tname, name, vol_frac) in zip(
                [:rock_thermal_conductivities, :fluid_thermal_conductivities],
                [:rock_thermal_conductivity, :fluid_thermal_conductivity],
                [1 .- ϕ, ϕ]
            )
            Λ = fractures[name, Cells()].*vol_frac
            T_hf = compute_half_face_trans(fmesh, 
            geo.cell_centroids, geo.face_centroids, 
            fractures[:areas, Faces()], Λ, faces, facepos)
            T = Jutul.EmbeddedMeshes.compute_face_trans_dfm(T_hf, N, fmesh.intersections)
            if contains(String(name), "fluid")
                if Λ isa Vector
                    T = repeat(T', number_of_phases(model[:Reservoir].system), 1)
                else
                    error("Fracture thermal conductivity $name must be a Vector.")
                end
            end
            fractures[tname, Faces()] = T
        end
    end

    fmodel = setup_reservoir_model(fractures, system; context = model.context, block_backend = block_backend, kwarg...)
    if fmodel isa Jutul.MultiModel
        fmodel = fmodel.models[:Reservoir]
    end
    fmc = JutulDarcy.FractureMatrixConnection()
    fmodel.domain.entities[fmc] = count_entities(fractures, fmc)

    # Reconstruct multimodel with fracture model and same non-fracture models as before
    new_models = JutulStorage()
    new_models[:Reservoir] = model.models[:Reservoir]
    new_models[:Fractures] = fmodel
    groups = Int[1, 1]
    for (k, v) in pairs(model.models)
        if k != :Reservoir
            new_models[k] = v
            push!(groups, 2)
        end
    end
    if isnothing(model.groups)
        groups = nothing
    end
    old_model = model
    model = Jutul.MultiModel(new_models; groups = groups)
    for ct in old_model.cross_terms
        push!(model.cross_terms, ct)
    end

    # Set up DFM cross-terms
    set_parameters!(model.models[:Fractures],
        FractureMatrixTransmissibility = FractureMatrixTransmissibility(),
        FractureMatrixGravityDifference = FractureMatrixGravityDifference(),
    )
    matrix_cells = fractures[:matrix_cells, fmc]
    fracture_cells = fractures[:connection_cells, fmc]
    ct = MatrixFromFractureFlowCT(matrix_cells, fracture_cells)
    add_cross_term!(model, ct, target = :Reservoir, source = :Fractures, equation = :mass_conservation)

    if has_thermal
        set_parameters!(model.models[:Fractures],
            FractureMatrixThermalConductivity = FractureMatrixThermalConductivity(),
        )
        ct = MatrixFromFractureThermalCT(matrix_cells, fracture_cells)
        add_cross_term!(model, ct, target = :Reservoir, source = :Fractures, equation = :energy_conservation)
    end

    for (name, well_model) in get_model_wells(model)
        g = physical_representation(well_model)
        if g isa WellDomain
            # WI = vec(g.perforations.WI)
            if !haskey(g.perforations, :fracture)
                continue
            end
            set_parameters!(model.models[name],
                FractureWellIndices = FractureWellIndices()
            )
            fc = vec(g.perforations.fracture)
            wc = vec(g.perforations.self_fracture)
            ct = FracturesFromWellFlowCT(fc, wc)
            add_cross_term!(model, ct, target = :Fractures, source = name, equation = :mass_conservation)
            if has_thermal
                ct = FracturesFromWellThermalCT(fc, wc)
                add_cross_term!(model, ct, target = :Fractures, source = name, equation = :energy_conservation)
                set_parameters!(model.models[name],
                    FractureWellIndicesThermal = FractureWellIndicesThermal()
                )
            end
        end
    end
    
    return model

end