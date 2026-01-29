using LinearAlgebra
using StaticArrays

# # Fix for Jutul.LinearizedBlock initialization with SMatrix
# function Base.convert(::Type{SMatrix{N, M, T, L}}, x::Number) where {N, M, T, L}
#     return SMatrix{N, M, T, L}(ntuple(_ -> T(x), Val(L)))
# end

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

    domain = reservoir_domain(mesh,
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

    domain[:matrix_faces] = matrix_faces

    return domain

end

function setup_reservoir_fracture_cross_term(reservoir::Jutul.DataDomain, fractures::Jutul.DataDomain)
    if ismissing(fractures[:matrix_faces])
        error("Fracture domain must have :matrix_faces data (from create_fracture_domain) to set up cross terms.")
    end
    matrix_faces = fractures[:matrix_faces]
    n_f = number_of_cells(fractures)
    n_res = number_of_cells(reservoir)
    
    # Prepare output arrays
    target_cells = Int64[]
    source_cells = Int64[]
    transmissibilities = Float64[]
    fmesh = physical_representation(fractures)

    vol_frac = fractures[:volumes]
    aperture = fractures[:aperture]

    mmesh = physical_representation(reservoir)
    N = get_neighborship(mmesh)

    for fcell in 1:n_f
        # Get area of fracture cell (face area)
        # Volume = Area * Aperture -> Area = Volume / Aperture
        a = fractures[:aperture][fcell]
        x_f = fractures[:cell_centroids][:, fcell]
        A_f = vol_frac[fcell]/a
        n_f = Jutul.EmbeddedMeshes.cell_normal(fmesh, fcell)
        K_f = fractures[:permeability][fcell]
        c_f = n_f.*aperture[fcell]/2.0
        T_fm = Jutul.half_face_trans(A_f, K_f, c_f, n_f)
        for face in matrix_faces[fcell]
            for cell in N[:, face]
                A_m = A_f
                n_m = n_f
                x_m = reservoir[:cell_centroids][:, cell]
                c_m = x_f .- x_m
                if dot(c_m, n_m) < 0
                    n_m = .-n_m
                end
                K_m = reservoir[:permeability][cell]
                T_mf = Jutul.half_face_trans(A_m, K_m, c_m, n_m)
                println("T_fm: $T_fm, T_mf: $T_mf")
                T = 1.0/(1.0/T_fm + 1.0/T_mf)
                push!(target_cells, cell)
                push!(source_cells, fcell)
                push!(transmissibilities, T)
            end
        end
    end
    
    return MatrixFromFractureFlowCT(target_cells, source_cells, transmissibilities)
end

function setup_reservoir_fracture_cross_term(reservoir::Jutul.SimulationModel, fractures::Jutul.SimulationModel)
    return setup_reservoir_fracture_cross_term(reservoir.data_domain, fractures.data_domain)
end

function JutulDarcy.setup_reservoir_model(reservoir::DataDomain, fractures::DataDomain, system::JutulSystem;
    wells = [],
    block_backend = true,
    kwarg...)

    # Set zero transmissibility accross matrix cells that are fractures
    T = reservoir_transmissibility(reservoir)
    T[fractures[:matrix_faces]] .= 0.0
    reservoir[:transmissibilities, Faces()] = T

    # Adjust matrix cell volumes to account for fracture volumes
    fracture_volumes = fractures[:volumes, Cells()]
    matrix_volumes = reservoir[:volumes, Cells()]
    N = get_neighborship(physical_representation(reservoir))
    for (cf, face) in enumerate(fractures[:matrix_faces])
        vf = fracture_volumes[cf]
        matrix_volumes[N[:, face]] .-= vf/2
    end
    reservoir[:volumes, Cells()] = matrix_volumes
    
    model = setup_reservoir_model(reservoir, system; wells = wells, block_backend = block_backend, kwarg...)
    fmodel = setup_reservoir_model(fractures, system; context = model.context, block_backend = false, kwarg...)
    if fmodel isa Jutul.MultiModel
        fmodel = fmodel.models[:Reservoir]
    end
    
    model.models[:Fractures] = fmodel
    if !isnothing(model.groups)
        group = maximum(model.groups)# + 1 # TODO: Now it gets schur lumped with the wells...
        push!(model.groups, group)
        model.group_lookup[:Fractures] = group
    end
    # Set up DFM cross-terms
    ct = setup_reservoir_fracture_cross_term(reservoir, fractures)
    add_cross_term!(model, ct, target = :Reservoir, source = :Fractures, equation = :mass_conservation)

    return model

end