using LinearAlgebra
using StaticArrays

function fracture_domain(mesh::Jutul.EmbeddedMeshes.EmbeddedMesh;
    aperture = 0.5e-3si_unit(:meter),
    permeability = missing,
    porosity = 1.0,
    kwarg...
    )

    all(isfinite, aperture) || throw(ArgumentError("Keyword argument aperture has non-finite entries."))
    minimum(aperture) >= 0 || throw(ArgumentError("All aperture values must be non-negative."))

    if ismissing(permeability)
        permeability = (aperture.^2)./12.0
    end

    domain = reservoir_domain(mesh;
        permeability = permeability,
        porosity = porosity,
        kwarg...
    )
    domain[:aperture, Cells()] = aperture

    geo = tpfv_geometry(mesh)
    aperture = domain[:aperture, Cells()]
    volumes = geo.volumes.*aperture
    volumes[mesh.intersection_cells] .*= aperture[mesh.intersection_cells]
    domain[:volumes, Cells()] = volumes

    N = domain[:neighbors]
    face_aperture = vec(sum(domain[:aperture, Cells()][N], dims = 1)./2)
    domain[:areas, Faces()] = geo.areas.*face_aperture
 
    star_delta = isempty(mesh.intersection_cells)
    nc = number_of_cells(mesh)
    faces, facepos = get_facepos(N, nc)
    T_hf = compute_half_face_trans(mesh, 
        geo.cell_centroids, geo.face_centroids, 
        domain[:areas, Faces()], permeability, aperture, faces, facepos)
    T = Jutul.EmbeddedMeshes.compute_face_trans_dfm(
            T_hf, N, mesh.intersection_neighbors, star_delta)
    domain[:transmissibilities, Faces()] = T

    return domain

end

function setup_matrix_fracture_cross_term(matrix::Jutul.DataDomain, fractures::Jutul.DataDomain, field::Symbol = :permeability)

    if field == :permeability
        matrix_conductivity = matrix[:permeability]
        fracture_conductivity = fractures[:permeability]
    elseif field == :thermal_conductivity
        function effective_conductivity(ϕ, Λ_r, Λ_f)
            return ϕ.*Λ_f .+ (1 .- ϕ).*Λ_r
        end
        matrix_conductivity = effective_conductivity(
            matrix[:porosity],
            matrix[:rock_thermal_conductivity],
            matrix[:fluid_thermal_conductivity]
        )
        fracture_conductivity = effective_conductivity(
            fractures[:porosity],
            fractures[:rock_thermal_conductivity],
            fractures[:fluid_thermal_conductivity]
        )
    else
        error("Unsupported field $field for matrix-fracture cross term setup.")
    end

    return setup_matrix_fracture_cross_term(matrix, fractures, matrix_conductivity, fracture_conductivity)
end

function setup_matrix_fracture_cross_term(matrix::Jutul.DataDomain, fractures::Jutul.DataDomain, matrix_conductivity, fracture_conductivity)
    if ismissing(fractures[:matrix_faces])
        error("Fracture domain must have :matrix_faces data (from create_fracture_domain) to set up cross terms.")
    end
    matrix_faces = fractures[:matrix_faces, FractureMatrixConnection()]
    n_frac = number_of_cells(fractures)
    n_res = number_of_cells(matrix)
    
    # Prepare output arrays
    target_cells = Int64[]
    source_cells = Int64[]
    transmissibilities = Float64[]
    gdz = Float64[]
    g = gravity_constant
    fmesh = physical_representation(fractures)

    vol_frac = fractures[:volumes]
    aperture = fractures[:aperture]

    mmesh = physical_representation(matrix)
    N = get_neighborship(mmesh)

    for fcell in 1:n_frac
        if fcell ∈ fmesh.intersection_cells
            continue # Skip intersection cells, as they are handled separately
        end
        # Get area of fracture cell (face area)
        # Volume = Area * Aperture -> Area = Volume / Aperture
        a = fractures[:aperture][fcell]
        x_f = fractures[:cell_centroids][:, fcell]
        A_f = vol_frac[fcell]/a
        n_f = Jutul.EmbeddedMeshes.cell_normal(fmesh, fcell)
        K_f = fracture_conductivity[fcell]
        c_f = n_f.*aperture[fcell]/2.0
        T_fm = Jutul.half_face_trans(A_f, K_f, c_f, n_f)
        for face in matrix_faces[fcell]
            for cell in N[:, face]
                A_m = A_f
                n_m = n_f
                x_m = matrix[:cell_centroids][:, cell]
                c_m = x_f .- x_m
                if dot(c_m, n_m) < 0
                    n_m = .-n_m
                end
                K_m = matrix_conductivity[cell]
                T_mf = Jutul.half_face_trans(A_m, K_m, c_m, n_m)
                T = 1.0/(1.0/T_fm + 1.0/T_mf)
                push!(target_cells, cell)
                push!(source_cells, fcell)
                push!(transmissibilities, T)
                dz = x_m[3] - x_f[3]
                push!(gdz, g*dz)
            end
        end
    end
    
    return target_cells, source_cells, transmissibilities, gdz
    
end

function compute_connection_transmissibilities(data_domain, source_cells, matrix_conductivity, matrix_centroids, fracture_conductivity)
    n_conn = length(source_cells)
    fmesh = physical_representation(data_domain)
    vol_frac = data_domain[:volumes]
    aperture = data_domain[:aperture]

    T = zeros(n_conn)
    for i in 1:n_conn
        fcell = source_cells[i]
        a = aperture[fcell]
        x_f = data_domain[:cell_centroids][:, fcell]
        A_f = vol_frac[fcell]/a
        n_f = Jutul.EmbeddedMeshes.cell_normal(fmesh, fcell)
        K_f = fracture_conductivity[fcell]
        c_f = n_f .* a/2.0
        T_fm = Jutul.half_face_trans(A_f, K_f, c_f, n_f)

        x_m = matrix_centroids[:, i]
        c_m = x_f .- x_m
        n_m = n_f
        if dot(c_m, n_m) < 0
            n_m = .-n_m
        end
        K_m = matrix_conductivity[i]
        T_mf = Jutul.half_face_trans(A_f, K_m, c_m, n_m)
        T[i] = 1.0/(1.0/T_fm + 1.0/T_mf)
    end
    return T
end

function JutulDarcy.setup_reservoir_model(matrix::DataDomain, fractures::DataDomain, system::Union{JutulSystem, Symbol};
    matrix_faces=missing,
    matrix_cells=missing,
    wells = [],
    block_backend = true,
    kwarg...)

    block_lump_with_res = true && block_backend

    mmesh = physical_representation(matrix)
    fmesh = physical_representation(fractures)
    if ismissing(matrix_faces)
        if haskey(fmesh, :parent_faces)
            matrix_faces = fmesh.parent_faces
        else
            @warn "Matrix faces not provided and not found in fracture mesh.\
            Fractures will not block affect through matrix cells."
        end
    end
    if ismissing(matrix_cells)
        !ismissing(matrix_faces) || error("Must provide matrix_cells if \
            matrix_faces is not provided to set up matrix/fracture connections.")
        N = get_neighborship(mmesh)
        matrix_cells = [c for c in eachcol(N[:, matrix_faces])]
    end
    fmc = FractureMatrixConnection()
    fmc = JutulDarcy.FractureMatrixConnection()
    num_fmc = length(matrix_faces)
    domain.entities[fmc] = num_fmc

    domain[:matrix_faces, fmc] = matrix_faces

    matrix_faces = fractures[:matrix_faces, fmc]
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
    
    model = setup_reservoir_model(matrix, system; wells = wells, block_backend = block_backend, kwarg...)
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

    fmodel = setup_reservoir_model(fractures, system; context = model.context, block_backend = block_backend && block_lump_with_res, kwarg...)
    if fmodel isa Jutul.MultiModel
        fmodel = fmodel.models[:Reservoir]
    end

    if block_lump_with_res
        if true
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
        else
            model.models[:Fractures] = fmodel
            if !isnothing(model.groups)
                push!(model.groups, 0)
                println(model.groups)
                for (k, (name, _)) in enumerate(pairs(model.models))
                    if name ∈ [:Reservoir, :Fractures]
                        gno = 1
                    else
                        gno = 2
                    end
                    model.groups[k] = gno
                    model.group_lookup[name] = gno
                end
            end
        end
    else
        model.models[:Fractures] = fmodel
        if !isnothing(model.groups)
            group = maximum(model.groups)# + 1 # TODO: Now it gets schur lumped with the wells...
            push!(model.groups, group)
            model.group_lookup[:Fractures] = group
        end
    end

    # Set up DFM cross-terms
    set_parameters!(model.models[:Fractures],
        FractureMatrixTransmissibility = FractureMatrixTransmissibility()
    )
    target_cells, source_cells, transmissibilities, gdz = setup_matrix_fracture_cross_term(matrix, fractures, :permeability)
    ct = MatrixFromFractureFlowCT(target_cells, source_cells, transmissibilities, gdz)
    add_cross_term!(model, ct, target = :Reservoir, source = :Fractures, equation = :mass_conservation)

    # Store per-connection matrix parameters on FractureMatrixConnection entities
    # fractures[:fracture_cells, fmc] = source_cells
    # fractures[:matrix_cell_centroids, fmc] = matrix[:cell_centroids][:, target_cells]
    # fractures[:matrix_permeability, fmc] = matrix[:permeability][target_cells]
    # fractures[:matrix_porosity, fmc] = matrix[:porosity][target_cells]

    if has_thermal
        target_cells, source_cells, transmissibilities_th, gdz = setup_matrix_fracture_cross_term(matrix, fractures, :thermal_conductivity)
        ct = MatrixFromFractureThermalCT(target_cells, source_cells, transmissibilities, transmissibilities_th, gdz)
        add_cross_term!(model, ct, target = :Reservoir, source = :Fractures, equation = :energy_conservation)
        fractures[:matrix_rock_thermal_conductivity, fmc] = matrix[:rock_thermal_conductivity][target_cells]
        fractures[:matrix_fluid_thermal_conductivity, fmc] = matrix[:fluid_thermal_conductivity][target_cells]
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