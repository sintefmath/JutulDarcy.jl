export compute_half_face_trans, compute_face_trans

function compute_half_face_trans(g, perm)
    compute_half_face_trans(g.cell_centroids, g.face_centroids, g.normals, g.areas, perm, g.neighbors)
end

function compute_half_face_trans(cell_centroids, face_centroids, face_normals, face_areas, perm, N)
    nf = size(N, 2)
    nhf = 2*nf
    dim = size(cell_centroids, 1)

    T_hf = similar(cell_centroids, nhf)
    faces, facePos = get_facepos(N)
    nc = length(facePos)-1
    if isa(perm, AbstractFloat)
        perm = repeat([perm], 1, nc)
    else
        perm::AbstractVecOrMat
    end

    # Sanity check
    @assert(dim == 2 || dim == 3)
    # Check cell centroids
    @assert(size(cell_centroids, 1) == dim)
    @assert(size(cell_centroids, 2) == nc)
    # Check face centroids
    @assert(size(face_centroids, 1) == dim)
    @assert(size(face_centroids, 2) == nf)
    # Check normals
    @assert(size(face_normals, 1) == dim)
    @assert(size(face_normals, 2) == nf)
    # Check areas
    @assert(length(face_areas) == nf)
    # Check perm
    @assert(size(perm, 2) == nc)
    # Check N, just in case
    @assert(size(N, 2) == nf)
    Threads.@threads for cell = 1:nc
        @inbounds for fpos = facePos[cell]:(facePos[cell+1]-1)
            face = faces[fpos]

            A = face_areas[face]
            K = expand_perm(perm[:, cell], dim)
            C = face_centroids[:, face] - cell_centroids[:, cell]
            Nn = face_normals[:, face]
            if N[2, face] == cell
                Nn = -Nn
            end
            T_hf[fpos] = compute_half_face_trans(A, K, C, Nn)
        end
    end
    return T_hf
end

function expand_perm(K, dim)
    if length(K) == dim
        return diagm(K)
    else
        @assert(length(K) == 1)
        return K[1]
    end
end

function compute_half_face_trans(A, K, C, N)
    return A*(dot(K*C, N))/dot(C, C)
end

function compute_face_trans(T_hf, N)
    faces, facePos = get_facepos(N)
    nf = size(N, 2)
    T = zeros(nf)
    for i in eachindex(faces)
        T[faces[i]] += 1/T_hf[i]
    end
    T = 1 ./T
    return T
    # (N[:, 2] .> 0) .& (N[:, 1] .> 0)
end

function discretized_domain_tpfv_flow(geometry; porosity = 0.1, 
                                                permeability = 9.869232667160131e-14, # 100 mD 
                                                T = nothing,
                                                fuse_flux = false,
                                                gravity = true,
                                                pore_volume = nothing)
    N = geometry.neighbors
    if isnothing(pore_volume)
        pore_volume = porosity.*geometry.volumes
    end
    if isnothing(T)
        T_hf = compute_half_face_trans(geometry, permeability)
        T = compute_face_trans(T_hf, N)
    end

    G = MinimalTPFAGrid(pore_volume, N)
    if dim(geometry) == 3 && gravity
        z = geometry.cell_centroids[3, :]
        g = gravity_constant
    else
        z = nothing
        g = nothing
    end

    if fuse_flux
        ft = DarcyMassMobilityFlowFused()
    else
        ft = DarcyMassMobilityFlow()
    end
    flow = TwoPointPotentialFlow(SPU(), TPFA(), ft, G, T, z, g, ncells = length(pore_volume))
    disc = (mass_flow = flow,)
    return DiscretizedDomain(G, disc)
end