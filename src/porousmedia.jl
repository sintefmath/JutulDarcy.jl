import Jutul: compute_half_face_trans, compute_face_trans

export compute_peaceman_index
function compute_peaceman_index(g::T, K, r, pos; kwarg...) where T<:Jutul.JutulMesh
    Δ = Jutul.cell_dims(g, pos)
    K = Jutul.expand_perm(K, dim(g), full = true)
    return compute_peaceman_index(Δ, K, r; kwarg...)
end

function compute_peaceman_index(Δ, K, radius; dir::Symbol = :z, constant = 0.14, Kh = nothing, skin = 0, check = true)
    K_d = diag(K)
    if dir == :x || dir == :X
        L, d1, d2 = Δ
        i, j = 2, 3
    elseif dir == :y || dir == :Y
        d1, L, d2 = Δ
        i, j = 1, 3
    else
        d1, d2, L = Δ
        i, j = 1, 2
        @assert dir == :z || dir == :Z "dir must be either :x, :y or :z (was :$dir)"
    end
    @assert L > 0
    @assert d1 > 0
    @assert d2 > 0
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

    if isnothing(Kh) || isnan(Kh)
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

function Jutul.discretize_domain(d::DataDomain, system::MultiPhaseSystem, ::Val{:default}; kwarg...)
    return discretized_domain_tpfv_flow(d; kwarg...)
end

export discretized_domain_tpfv_flow
function discretized_domain_tpfv_flow(geometry; porosity = 0.1, 
                                                permeability = 9.869232667160131e-14, # 100 mD 
                                                T = nothing,
                                                gdz = nothing,
                                                gravity = true,
                                                pore_volume = nothing,
                                                general_ad = false)
    N = geometry.neighbors
    error("Deprecated.")
    if isnothing(pore_volume)
        pore_volume = porosity.*geometry.volumes
    end
    nc = length(pore_volume)
    if isnothing(T)
        T = compute_face_trans(geometry, permeability)
    end
    cc = geometry.cell_centroids
    if gravity && size(cc, 1) == 3
        z = vec(cc[3, :])
    else
        z = zeros(nc)
    end
    if isnothing(gdz)
        gdz = compute_face_gdz(N, z)
    end
    G = MinimalTPFAGrid(pore_volume, N, trans = T, gdz = gdz)
    if general_ad
        flow = PotentialFlow(N, nc)
    else
        flow = TwoPointPotentialFlowHardCoded(G, ncells = nc)
    end
    disc = (mass_flow = flow, heat_flow = flow)
    return DiscretizedDomain(G, disc)
end

function discretized_domain_tpfv_flow(domain::Jutul.DataDomain; general_ad = false)
    N = domain[:neighbors]
    nc = number_of_cells(physical_representation(domain))
    if general_ad
        d = PotentialFlow(N, nc)
    else
        d = TwoPointPotentialFlowHardCoded(N, nc)
    end
    disc = (mass_flow = d, heat_flow = d)
    G = MinimalTPFATopology(N, ncells = nc)
    return DiscretizedDomain(G, disc)
end

function Jutul.discretize_domain(d::DataDomain{W}, system::MultiPhaseSystem, ::Val{:default}; kwarg...) where W<:Union{SimpleWell, MultiSegmentWell}
    return discretized_domain_well(physical_representation(d); kwarg...)
end

export discretized_domain_well
function discretized_domain_well(W::MultiSegmentWell; z = nothing, kwarg...)
    if isnothing(z)
        z = vec(W.centers[3, :])
    end
    flow = WellSegmentFlow(W, z)
    disc = (mass_flow = flow, heat_flow = flow)
    return DiscretizedDomain(W, disc; kwarg...)
end

function discretized_domain_well(W::SimpleWell; z = nothing, kwarg...)
    disc = (mass_flow = PotentialFlow(W), heat_flow = PotentialFlow(W))
    return DiscretizedDomain(W, disc; kwarg...)
end
