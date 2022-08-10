import Jutul: compute_half_face_trans, compute_face_trans

export compute_peaceman_index
function compute_peaceman_index(g::T, K, r, pos; kwarg...) where T<:Jutul.AbstractJutulMesh
    Δ = Jutul.cell_dims(g, pos)
    K = Jutul.expand_perm(K, dim(g), full = true)
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
                                                gravity = true,
                                                pore_volume = nothing)
    N = geometry.neighbors
    if isnothing(pore_volume)
        pore_volume = porosity.*geometry.volumes
    end
    if isnothing(T)
        T = compute_face_trans(geometry, permeability)
    end

    G = MinimalTPFAGrid(pore_volume, N)
    if dim(geometry) == 3 && gravity
        z = geometry.cell_centroids[3, :]
        g = gravity_constant
    else
        z = nothing
        g = nothing
    end

    flow = TwoPointPotentialFlowHardCoded(G, T, z, g, ncells = length(pore_volume))
    disc = (mass_flow = flow,)
    return DiscretizedDomain(G, disc)
end

export discretized_domain_well
function discretized_domain_well(W; z = nothing, kwarg...)
    if W isa MultiSegmentWell
        if isnothing(z)
            z = vec(W.centers[3, :])
        end
    else
        z = []
    end
    flow = WellSegmentFlow(W, z)
    disc = (mass_flow = flow,)
    return DiscretizedDomain(W, disc; kwarg...)
end
