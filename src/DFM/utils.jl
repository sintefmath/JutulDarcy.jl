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

    return domain

end

function JutulDarcy.setup_reservoir_model(reservoir::DataDomain, fractures::DataDomain, system::JutulSystem;
    wells = [],
    kwarg...)

    model = setup_reservoir_model(reservoir, system; wells = wells, kwarg...)
    fmodel = setup_reservoir_model(fractures, system; kwarg...)
    model.models[:Fractures] = fmodel

    # Set up DFM cross-terms

    return model

end