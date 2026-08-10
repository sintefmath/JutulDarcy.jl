# function add_thermal_dispersion(domain::DataDomain, D::AbstractVector{<:Real})x
    
#     D = if length(D) == 2
#         repeat(D, 1, number_of_cells(domain))
#     elseif length(D) == 1
#         domain[:thermal_dispersivity] = repeat([D[1], D[1]], 1, number_of_cells(domain))
#     else
#         error("Unsupported size for thermal dispersivity array")
#     end

# end


function add_thermal_dispersion!(domain::DataDomain, D::AbstractArray{<:Real, 2})

    # Add to data domain based on size of D
    D = if size(D, 1) == 2 && size(D, 2) == 1
        repeat(D, 1, number_of_cells(domain))
    elseif size(D, 1) == 1 && size(D, 2) == number_of_cells(domain)
        repeat(D, 2, 1)
    else
        error("Unsupported size for thermal dispersivity array")
    end
    domain[:thermal_dispersivity] = D

    # Add useful geometric quantities
    neighbors = domain[:neighbors]
    nc = size(domain[:cell_centroids], 2)
    faces, facepos = get_facepos(neighbors, nc)
    facesigns = Jutul.get_facesigns(neighbors, faces, facepos, nc)

    half_face_vectors = zeros(size(domain[:cell_centroids], 1), length(faces))

    for cell in 1:number_of_cells(domain)
        xc = domain[:cell_centroids][:, cell]
        for fpos = facepos[cell]:(facepos[cell+1]-1)
            face = faces[fpos]
            sgn = facesigns[fpos]
            xf = domain[:face_centroids][:, face]
            half_face_vectors[:, fpos] = (xf .- xc)
        end
    end
    hf = HalfFaces()
    domain[:half_faces, hf] = faces
    domain[:half_facepos, NoEntity()] = facepos
    domain[:half_facesigns, hf] = facesigns
    domain[:half_face_vectors, hf] = half_face_vectors
    
    return domain

end