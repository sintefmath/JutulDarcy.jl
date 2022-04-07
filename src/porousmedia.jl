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

function expand_perm(K, dim; full = false)
    n = length(K)
    if n == dim
        K_e = diagm(K)
    elseif n == 1
        K_e = first(K)
        if full
            # Expand to matrix
            K_e = zeros(dim, dim) + I*K_e
        end
    else
        if dim == 2
            @assert n == 3 "Two-dimensional grids require 1/2/3 permeability entries per cell (was $n)"
            K_e = [K[1] K[2]; K[2] K[3]]
        else
            @assert n == 6 "Three-dimensional grids require 1/3/6 permeability entries per cell (was $n)"
            K_e =  [K[1] K[2] K[3];
                    K[2] K[4] K[5];
                    K[3] K[5] K[6]]
        end
    end
    return K_e
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
end

export compute_peaceman_index
function compute_peaceman_index(g::T, K, r, pos; kwarg...) where T<:Jutul.AbstractJutulMesh
    Δ = Jutul.cell_dims(g, pos)
    K = expand_perm(K, dim(g), full = true)
    return compute_peaceman_index(Δ, K, r; kwarg...)
end

function compute_peaceman_index(Δ, K, radius; dir::Symbol = :z, constant = 0.14, Kh = nothing, skin = 0, check = true)
    K_d = diag(K)
    if dir == :x
        L, d1, d2 = Δ
        i, j = 2, 3
    elseif dir == :y
        d1, L, d2 = Δ
        i, j = 1, 3
    else
        d1, d2, L = Δ
        i, j = 1, 2
        @assert dir == :z "dir must be either :x, :y or :z (was :$dir)"
    end
    k1, k2 = K_d[i], K_d[j]

    function kratio(l, v)
        r = l/v
        if isfinite(r)
            return r
        else
            return zero(l)
        end
    end
    k21 = kratio(k2, k1)
    k12 = kratio(k1, k2)

    re1 = 2 * constant * sqrt((d1^2)*sqrt(k21) + (d2^2)*sqrt(k12))
    re2 = k21^(1/4) + k12^(1/4)
 
    re  = kratio(re1, re2)
    ke  = sqrt(k1*k2)
 
    if isnothing(Kh)
        Kh = L*ke
    end
    WI = 2 * π * Kh / (log(re / radius) + skin)
    if check && WI < 0
        if re < radius
            error("Equivialent Peaceman radius is smaller than well radius - negative well was negative. Either the cell is too small, or the radius too big.")
        else
            error("Too large skin factor - well radius became negative.")
        end
    end
    return WI
end

export discretized_domain_tpfv_flow
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