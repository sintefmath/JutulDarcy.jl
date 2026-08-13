
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
    domain[:face_fluxes_tmp, Faces()] = zeros(number_of_faces(domain))

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


Base.@kwdef mutable struct TotalDarcyVelocity <: VectorVariables
    change_tol_rel = 1e-3
    changed = true
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
    nf = number_of_faces(domain)
    volumes = domain[:volumes]

    half_facesigns = domain[:half_facesigns]
    half_face_vectors_svec = domain[:half_face_vectors_svec]
    facepos = domain[:half_facepos]
    faces = domain[:half_faces]

    face_fluxes = domain[:face_fluxes_tmp]
    neighbors = domain[:neighbors_svec]
    @inbounds for face in 1:nf
        face_fluxes[face] = face_total_volume_flux(face, neighbors, state, model)
    end

    v0 = copy(v)
    change = false
    @views fill!(v, 0.0)
    @inbounds for cell in 1:nc
        v_c = v[:, cell]
        for fpos = facepos[cell]:(facepos[cell + 1] - 1)
            face = faces[fpos]
            q = face_fluxes[face]*half_facesigns[fpos]
            if q != 0.0
                v_c .+= q.*half_face_vectors_svec[fpos]
            end
        end
        v_c ./= volumes[cell]
        change = change ||
            any(norm(v_c .- v0[:, cell]) .> prm.change_tol_rel .* norm(v0[:, cell]))
        v[:, cell] = v_c
    end
    prm.changed = change
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
        
function thermal_dispersion_half_face_transmissibility(face, side, cell_centroids, face_centroids, normals, neighbors, areas, v_cell, aT_cell, aL_cell, is_active_cell)
    @inbounds left, right = neighbors[face]
    if side <= 0 || right <= 0
        return 0.0
    end
    if !is_active_cell[side]
        return 0.0
    end
    @inbounds A = areas[face]
    @inbounds Nf = normals[face]
    @inbounds xf = face_centroids[face]

    v = v_cell[side]
    aT = aT_cell[side]
    aL = aL_cell[side]
    C = xf - cell_centroids[side]
    sgn = ifelse(side == left, 1.0, -1.0)
    N = sgn*Nf
    # D*C where D = aT*I + aL*(v*v')
    KC = aT*C + aL*v*dot(v, C)
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

    if !model.parameters[:TotalDarcyVelocity].changed
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
    nc = number_of_cells(domain)

    # Precompute cell-wise velocity-dependent dispersion coefficients.
    v_cell = Vector{eltype(cell_centroids)}(undef, nc)
    aT_cell = zeros(nc)
    aL_cell = zeros(nc)
    is_active_cell = falses(nc)
    @inbounds for cell in 1:nc
        v_cell_i = eltype(cell_centroids)(v[:, cell])
        v_cell[cell] = v_cell_i
        vmag = sqrt(dot(v_cell_i, v_cell_i))
        if vmag > 0.0
            α_L = state.ThermalDispersivity[1, cell]
            α_T = state.ThermalDispersivity[2, cell]
            aT_cell[cell] = α_T*vmag
            aL_cell[cell] = (α_L - α_T)/vmag
            is_active_cell[cell] = true
        end
    end

    @inbounds for face in 1:nf
        left, right = neighbors[face]
        θ_d[1, face] = thermal_dispersion_half_face_transmissibility(face, left, cell_centroids, face_centroids, normals, neighbors, areas, v_cell, aT_cell, aL_cell, is_active_cell)
        θ_d[2, face] = thermal_dispersion_half_face_transmissibility(face, right, cell_centroids, face_centroids, normals, neighbors, areas, v_cell, aT_cell, aL_cell, is_active_cell)
    end
    return θ_d
end
