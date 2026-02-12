using LinearAlgebra
using StaticArrays

function fracture_domain(mesh::Jutul.EmbeddedMeshes.EmbeddedMesh;
    aperture = 0.5e-3si_unit(:meter),
    permeability = missing,
    porosity = 1.0,
    matrix_faces = missing,
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
    domain[:volumes, Cells()] = geo.volumes.*domain[:aperture, Cells()]
    
    N = domain[:neighbors]
    face_aperture = vec(sum(domain[:aperture, Cells()][N], dims = 1)./2)
    domain[:areas, Faces()] = geo.areas.*face_aperture

    nc = number_of_cells(mesh)
    faces, facepos = get_facepos(N, nc)
    T_hf = compute_half_face_trans(mesh, 
        geo.cell_centroids, geo.face_centroids, 
        domain[:areas, Faces()], permeability, faces, facepos)
    T = Jutul.EmbeddedMeshes.compute_face_trans_dfm(T_hf, N, mesh.intersections)
    domain[:transmissibilities, Faces()] = T

    domain[:matrix_faces, NoEntity()] = matrix_faces

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
    matrix_faces = fractures[:matrix_faces, NoEntity()]
    n_f = number_of_cells(fractures)
    n_res = number_of_cells(matrix)
    
    # Prepare output arrays
    target_cells = Int64[]
    source_cells = Int64[]
    transmissibilities = Float64[]
    fmesh = physical_representation(fractures)

    vol_frac = fractures[:volumes]
    aperture = fractures[:aperture]

    mmesh = physical_representation(matrix)
    N = get_neighborship(mmesh)

    for fcell in 1:n_f
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
            end
        end
    end
    
    return target_cells, source_cells, transmissibilities
    
end

function JutulDarcy.setup_reservoir_model(matrix::DataDomain, fractures::DataDomain, system::Union{JutulSystem, Symbol};
    wells = [],
    block_backend = true,
    kwarg...)

    block_lump_with_res = true# && block_backend

    # Set zero transmissibility accross matrix cells that are fractures
    T = reservoir_transmissibility(matrix)
    T[fractures[:matrix_faces, NoEntity()]] .= 0.0
    matrix[:transmissibilities, Faces()] = T

    if haskey(matrix, :rock_thermal_conductivity)
        T = reservoir_conductivity(matrix)
        T[fractures[:matrix_faces, NoEntity()]] .= 0.0
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
        T[fractures[:matrix_faces, NoEntity()]] .= 0.0
        matrix[:fluid_thermal_conductivities, Faces()] = T
    end

    # Adjust matrix cell volumes to account for fracture volumes
    fracture_volumes = fractures[:volumes, Cells()]
    matrix_volumes = matrix[:volumes, Cells()]
    N = get_neighborship(physical_representation(matrix))
    for (cf, face) in enumerate(fractures[:matrix_faces, NoEntity()])
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
        model = Jutul.MultiModel(new_models; groups = groups)
    else
        model.models[:Fractures] = fmodel
        if !isnothing(model.groups)
            group = maximum(model.groups)# + 1 # TODO: Now it gets schur lumped with the wells...
            push!(model.groups, group)
            model.group_lookup[:Fractures] = group
        end
    end

    # Set up DFM cross-terms
    target_cells, source_cells, transmissibilities = setup_matrix_fracture_cross_term(matrix, fractures, :permeability)
    gdz = zeros(Float64, length(transmissibilities))
    ct = MatrixFromFractureFlowCT(target_cells, source_cells, transmissibilities, gdz)
    add_cross_term!(model, ct, target = :Reservoir, source = :Fractures, equation = :mass_conservation)
    if has_thermal
        target_cells, source_cells, transmissibilities_th = setup_matrix_fracture_cross_term(matrix, fractures, :thermal_conductivity)
        ct = MatrixFromFractureThermalCT(target_cells, source_cells, transmissibilities, transmissibilities_th, gdz)
        add_cross_term!(model, ct, target = :Reservoir, source = :Fractures, equation = :energy_conservation)
    end

    for (name, well_model) in get_model_wells(model)
        g = physical_representation(well_model)
        if g isa WellDomain
            # WI = vec(g.perforations.WI)
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