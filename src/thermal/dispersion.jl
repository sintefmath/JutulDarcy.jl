
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
            xf = domain[:face_centroids][:, face]
            half_face_vectors[:, fpos] = (xf .- xc)
        end
    end
    hf = HalfFaces()
    domain[:half_faces, hf] = faces
    domain[:half_facepos, NoEntity()] = facepos
    domain[:half_facesigns, hf] = facesigns
    domain[:half_face_vectors, hf] = half_face_vectors
    dim = size(domain[:cell_centroids], 1)
    domain[:half_face_vectors_svec, hf] = [SVector{dim}(half_face_vectors[:, i]) for i in axes(half_face_vectors, 2)]

    # Cache static geometry as SVectors for reuse in hot loops.
    domain[:cell_centroids_svec, Cells()] = [SVector{dim}(domain[:cell_centroids][:, i]) for i in axes(domain[:cell_centroids], 2)]
    domain[:face_centroids_svec, Faces()] = [SVector{dim}(domain[:face_centroids][:, i]) for i in axes(domain[:face_centroids], 2)]
    domain[:normals_svec, Faces()] = [SVector{dim}(domain[:normals][:, i]) for i in axes(domain[:normals], 2)]
    if neighbors isa AbstractMatrix
        domain[:neighbors_svec, Faces()] = [SVector{2, Int}(neighbors[1, i], neighbors[2, i]) for i in axes(neighbors, 2)]
    else
        domain[:neighbors_svec, Faces()] = [SVector{2, Int}(n[1], n[2]) for n in neighbors]
    end
    
    return domain

end


struct TotalDarcyVelocity <: VectorVariables end
# Jutul.variable_scale(::ThermalDispersion) = 1.0
# Jutul.minimum_value(::ThermalDispersion) = 0.0
Jutul.values_per_entity(model, ::TotalDarcyVelocity) = dim(model.data_domain)
Jutul.associated_entity(::TotalDarcyVelocity) = Cells()

function Jutul.default_parameter_values(data_domain, model, param::TotalDarcyVelocity, symb)
    dim = size(data_domain[:cell_centroids], 1)
    nc = number_of_cells(data_domain)
    return zeros(dim, nc)
end

function Jutul.update_parameter_before_step!(v, ::TotalDarcyVelocity, storage, model, dt, forces)
    state = storage.state
    domain = model.data_domain
    nc = number_of_cells(domain)
    nf = number_of_faces(domain)
    volumes = domain[:volumes]

    half_facesigns = domain[:half_facesigns]
    half_face_vectors_svec = domain[:half_face_vectors_svec]
    facepos = domain[:half_facepos]
    faces = domain[:half_faces]

    face_fluxes = zeros(nf)
    neighbors = domain[:neighbors_svec]
    @inbounds for face in 1:nf
        q_face = face_flux_helper(face, neighbors, state, model, false)
        q_sum = 0.0
        for q in q_face
            q_sum += value(q)
        end
        face_fluxes[face] = q_sum
    end

    @views fill!(v, 0.0)
    @inbounds for cell in 1:nc
        cell_velocity = v[:, cell]
        for fpos = facepos[cell]:(facepos[cell + 1] - 1)
            face = faces[fpos]
            q = face_fluxes[face]*half_facesigns[fpos]
            if q != 0.0
                cell_velocity .+= q.*half_face_vectors_svec[fpos]
            end
        end
        cell_velocity ./= volumes[cell]
        v[:, cell] = cell_velocity
    end
    return v
end

struct ThermalDispersionTransmissibility <: VectorVariables end
Jutul.minimum_value(::ThermalDispersionTransmissibility) = 0.0
Jutul.values_per_entity(model, ::ThermalDispersionTransmissibility) = 2
Jutul.associated_entity(::ThermalDispersionTransmissibility) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::ThermalDispersionTransmissibility, symb)
    return zeros(2, number_of_faces(data_domain))
end
        
function thermal_dispersion_half_face_transmissibility(face, side, state, cell_centroids, face_centroids, normals, neighbors, areas)
    @inbounds left, right = neighbors[face]
    if side <= 0 || right <= 0
        return 0.0
    end
    @inbounds A = areas[face]
    @inbounds Nf = normals[face]
    @inbounds xf = face_centroids[face]

    dim = size(cell_centroids[1], 1)
    Id = SMatrix{dim, dim, Float64}(I)
    if haskey(state, :TotalDarcyVelocity)
        v = @view state.TotalDarcyVelocity[:, side]
    else
        v = zeros(size(cell_centroids[1]))
    end
    v² = sqrt(v'*v)
    if v² == 0.0
        return 0.0
    end
    α_L = state.ThermalDispersivity[1, side]
    α_T = state.ThermalDispersivity[2, side]
    D = (α_T*Id + (α_L - α_T)*(v*v')/(v'*v)).*sqrt((v'*v))
    C = xf - cell_centroids[side]
    sgn = ifelse(side == left, 1.0, -1.0)
    return Jutul.half_face_trans(A, D, C, sgn*Nf)
end

function Jutul.update_parameter_before_step!(θ_d, ::ThermalDispersionTransmissibility, storage, model, dt, forces)
    state = storage.state
    if !haskey(state, :ThermalDispersivity)
        fill!(θ_d, 0.0)
        return θ_d
    end
    dispersion = state.ThermalDispersivity
    if all(x -> value(x) == 0.0, dispersion)
        fill!(θ_d, 0.0)
        return θ_d
    end

    domain = model.data_domain
    nf = number_of_faces(domain)
    cell_centroids = domain[:cell_centroids_svec]
    face_centroids = domain[:face_centroids_svec]
    normals = domain[:normals_svec]
    neighbors = domain[:neighbors_svec]
    areas = domain[:areas]

    @inbounds for face in 1:nf
        left, right = neighbors[face]
        θ_d[1, face] = thermal_dispersion_half_face_transmissibility(face, left, state, cell_centroids, face_centroids, normals, neighbors, areas)
        θ_d[2, face] = thermal_dispersion_half_face_transmissibility(face, right, state, cell_centroids, face_centroids, normals, neighbors, areas)
    end
    return θ_d
end
