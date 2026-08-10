
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
    half_face_vectors = domain[:half_face_vectors]
    facepos = domain[:half_facepos]
    faces = domain[:half_faces]

    face_fluxes = zeros(nf)
    neighbors = domain[:neighbors]
    @inbounds for face in 1:nf
        q_face = face_flux_helper(face, neighbors, state, model, false)
        face_fluxes[face] = sum(value.(q_face))
    end

    @views fill!(v, 0.0)
    @inbounds for cell in 1:nc
        cell_velocity = v[:, cell]
        for fpos = facepos[cell]:(facepos[cell + 1] - 1)
            face = faces[fpos]
            q = face_fluxes[face]*half_facesigns[fpos]
            if q != 0.0
                cell_velocity .+= q.*@view(half_face_vectors[:, fpos])
            end
        end
        cell_velocity ./= volumes[cell]
        v[:, cell] = cell_velocity
    end
    return v
end

struct ThermalDispersionTransmissibility <: ScalarVariable end
Jutul.minimum_value(::ThermalDispersionTransmissibility) = 0.0
Jutul.associated_entity(::ThermalDispersionTransmissibility) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::ThermalDispersionTransmissibility, symb)
    return zeros(number_of_faces(data_domain))
end

function to_vec_of_svectors(x::AbstractVector{<:SVector})
    return x
end

function to_vec_of_svectors(x::AbstractMatrix)
    dim = size(x, 1)
    return [SVector{dim}(x[:, i]) for i in axes(x, 2)]
end

function to_vec_of_svectors(x)
    return x
end

function thermal_dispersion_transmissibility(face, state, model)
    domain = model.data_domain
    cell_centroids = to_vec_of_svectors(domain[:cell_centroids])
    face_centroids = to_vec_of_svectors(domain[:face_centroids])
    normals = to_vec_of_svectors(domain[:normals])
    neighbors = to_vec_of_svectors(domain[:neighbors])
    areas = domain[:areas]
    component_heat_capacity = domain[:component_heat_capacity]
    return thermal_dispersion_transmissibility(face, state, model, cell_centroids, face_centroids, normals, neighbors, areas, component_heat_capacity)
end

function thermal_dispersion_transmissibility(face, state, model, cell_centroids, face_centroids, normals, neighbors, areas, component_heat_capacity)
    @inbounds left, right = neighbors[face]
    @inbounds A = areas[face]
    @inbounds Nf = normals[face]
    @inbounds xf = face_centroids[face]

    den = 0.0
    Id = SMatrix{3, 3, Float64}(I)
    for side in (left, right)
        if haskey(state, :TotalDarcyVelocity)
            v = @view state.TotalDarcyVelocity[:, side]
        else
            v = zeros(size(cell_centroids[1]))
        end
        v² = sqrt(v'*v)
        if v² == 0.0
            continue
        end
        α_L = state.ThermalDispersivity[1, side]
        α_T = state.ThermalDispersivity[2, side]
        ρ_f = value(state.PhaseMassDensities[1, side])
        Cp_f = component_heat_capacity[side]
        D = ρ_f*Cp_f*(α_T*Id + (α_L - α_T)*(v*v')/(v'*v)).*sqrt((v'*v))
        C = xf - cell_centroids[side]
        sgn = ifelse(side == left, 1.0, -1.0)
        θ_hf = Jutul.half_face_trans(A, D, C, sgn*Nf)
        den += 1/θ_hf
    end
    return den == 0.0 ? 0.0 : 1/den
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
    cell_centroids = to_vec_of_svectors(domain[:cell_centroids])
    face_centroids = to_vec_of_svectors(domain[:face_centroids])
    normals = to_vec_of_svectors(domain[:normals])
    neighbors = to_vec_of_svectors(domain[:neighbors])
    areas = domain[:areas]
    component_heat_capacity = domain[:component_heat_capacity]

    @inbounds for face in 1:nf
        θ_d[face] = thermal_dispersion_transmissibility(face, state, model, cell_centroids, face_centroids, normals, neighbors, areas, component_heat_capacity)
    end
    return θ_d
end
