
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

Base.@kwdef mutable struct TotalVolumeFlux <: ScalarVariable
    change_tol_rel = 1e-2
    changed = BitVector()
    changed_any = true
end
Jutul.associated_entity(::TotalVolumeFlux) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::TotalVolumeFlux, symb)
    return zeros(number_of_faces(data_domain))
end

function Jutul.update_parameter_before_step!(q_t, prm::TotalVolumeFlux, storage, model, dt, forces)
    state = storage.state
    domain = model.data_domain
    neighbors = domain[:neighbors_svec]
    nf = number_of_faces(domain)

    if length(prm.changed) != nf
        prm.changed = falses(nf)
    end
    changed_faces = prm.changed
    changed_any = false

    @inbounds for face in 1:nf
        q_old = q_t[face]
        q_new = face_total_volume_flux(face, neighbors, state, model)
        q_t[face] = q_new

        thresh = prm.change_tol_rel*max(abs(q_old), eps(Float64))
        has_changed = abs(q_new - q_old) > thresh
        changed_faces[face] = has_changed
        changed_any = changed_any || has_changed
    end
    prm.changed_any = changed_any
    return q_t
end


Base.@kwdef mutable struct TotalDarcyVelocity <: VectorVariables
    changed = BitVector()
    changed_any = true
end
# Jutul.variable_scale(::ThermalDispersion) = 1.0
# Jutul.minimum_value(::ThermalDispersion) = 0.0
Jutul.values_per_entity(model, ::TotalDarcyVelocity) = dim(model.data_domain)
Jutul.associated_entity(::TotalDarcyVelocity) = Cells()

function Jutul.default_parameter_values(data_domain, model, param::TotalDarcyVelocity, symb)
    dim = size(data_domain[:cell_centroids], 1)
    nc = number_of_cells(data_domain)
    return zeros(dim, nc)
end

function Jutul.update_parameter_before_step!(v, prm::TotalDarcyVelocity, storage, model, dt, forces)
    state = storage.state
    domain = model.data_domain
    nc = number_of_cells(domain)
    volumes = domain[:volumes]

    half_facesigns = domain[:half_facesigns]
    half_face_vectors_svec = domain[:half_face_vectors_svec]
    facepos = domain[:half_facepos]
    faces = domain[:half_faces]

    face_fluxes = state.TotalVolumeFlux

    if length(prm.changed) != nc
        prm.changed = falses(nc)
    end
    changed_cells = prm.changed
    fill!(changed_cells, false)

    flux_prm = model.parameters[:TotalVolumeFlux]
    changed_faces = flux_prm.changed
    if length(changed_faces) == 0
        fill!(changed_cells, true)
    else
        @inbounds for face in eachindex(changed_faces)
            if changed_faces[face]
                left, right = domain[:neighbors_svec][face]
                if left > 0
                    changed_cells[left] = true
                end
                if right > 0
                    changed_cells[right] = true
                end
            end
        end
    end

    changed_any = any(changed_cells)
    prm.changed_any = changed_any
    if !changed_any
        return v
    end

    @inbounds for cell in 1:nc
        if !changed_cells[cell]
            continue
        end
        v_c = v[:, cell]
        fill!(v_c, 0.0)
        for fpos = facepos[cell]:(facepos[cell + 1] - 1)
            face = faces[fpos]
            q = face_fluxes[face]*half_facesigns[fpos]
            if q != 0.0
                v_c .+= q.*half_face_vectors_svec[fpos]
            end
        end
        v_c ./= volumes[cell]
        v[:, cell] = v_c
    end
    return v
end

@inline function face_total_volume_flux(face, neighbors, state, model)
    q_face = face_flux_helper(face, neighbors, state, model, false)
    q_sum = 0.0
    @inbounds for q in q_face
        q_sum += value(q)
    end
    return q_sum
end

struct ThermalDispersionTransmissibility <: VectorVariables end
Jutul.minimum_value(::ThermalDispersionTransmissibility) = 0.0
Jutul.values_per_entity(model, ::ThermalDispersionTransmissibility) = 2
Jutul.associated_entity(::ThermalDispersionTransmissibility) = Faces()

function Jutul.default_parameter_values(data_domain, model, param::ThermalDispersionTransmissibility, symb)
    return zeros(2, number_of_faces(data_domain))
end
        
function thermal_dispersion_half_face_transmissibility(face, side, cell_centroids, face_centroids, normals, neighbors, areas, v, dispersion)
    @inbounds left, right = neighbors[face]
    if side <= 0 || right <= 0
        return 0.0
    end
    @inbounds A = areas[face]
    @inbounds Nf = normals[face]
    @inbounds xf = face_centroids[face]

    v_i = eltype(cell_centroids)(v[:, side])
    vmag = sqrt(dot(v_i, v_i))
    if vmag == 0.0
        return 0.0
    end
    α_L = value(dispersion[1, side])
    α_T = value(dispersion[2, side])
    aT = α_T*vmag
    aL = (α_L - α_T)/vmag

    C = xf - cell_centroids[side]
    sgn = ifelse(side == left, 1.0, -1.0)
    N = sgn*Nf
    # D*C where D = aT*I + aL*(v*v')
    KC = aT*C + aL*v_i*dot(v_i, C)
    return A*dot(KC, N)/dot(C, C)
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

    velocity_prm = model.parameters[:TotalDarcyVelocity]
    if !velocity_prm.changed_any
        return θ_d
    end

    domain = model.data_domain
    nf = number_of_faces(domain)
    cell_centroids = domain[:cell_centroids_svec]
    face_centroids = domain[:face_centroids_svec]
    normals = domain[:normals_svec]
    neighbors = domain[:neighbors_svec]
    areas = domain[:areas]
    v = state.TotalDarcyVelocity
    changed_cells = velocity_prm.changed

    @inbounds for face in 1:nf
        left, right = neighbors[face]
        left_changed = left > 0 && changed_cells[left]
        right_changed = right > 0 && changed_cells[right]
        if !(left_changed || right_changed)
            continue
        end
        θ_d[1, face] = thermal_dispersion_half_face_transmissibility(face, left, cell_centroids, face_centroids, normals, neighbors, areas, v, dispersion)
        θ_d[2, face] = thermal_dispersion_half_face_transmissibility(face, right, cell_centroids, face_centroids, normals, neighbors, areas, v, dispersion)
    end
    return θ_d
end
