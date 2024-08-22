import Jutul: compute_half_face_trans, compute_face_trans

function compute_peaceman_index(g::T, K, r, pos; kwarg...) where T<:Jutul.JutulMesh
    Δ = peaceman_cell_dims(g, pos)
    K = Jutul.expand_perm(K, dim(g))
    return compute_peaceman_index(Δ, K, r; kwarg...)
end

function peaceman_cell_dims(g, pos)
    horz = get_mesh_entity_tag(g, Faces(), :orientation, :horizontal, throw = false)
    vert = get_mesh_entity_tag(g, Faces(), :orientation, :vertical, throw = false)
    if !(g isa UnstructuredMesh) || ismissing(horz) || ismissing(vert) || Jutul.dim(g) < 3
        Δ = Jutul.cell_dims(g, pos)
    else
        index = cell_index(g, pos)
        xy_min = SVector{2, Float64}(Inf, Inf)
        xy_max = SVector{2, Float64}(-Inf, -Inf)

        z_min = Inf
        z_max = -Inf
        for (e, face_set) in [(Faces(), g.faces), (BoundaryFaces(), g.boundary_faces)]
            for face in face_set.cells_to_faces[index]
                face_centroid, = Jutul.compute_centroid_and_measure(g, e, face)
                if mesh_entity_has_tag(g, e, :orientation, :horizontal, face)
                    z_min = min.(z_min, face_centroid[3])
                    z_max = max.(z_max, face_centroid[3])
                elseif mesh_entity_has_tag(g, e, :orientation, :vertical, face)
                    xy_min = min.(xy_min, face_centroid[1:2])
                    xy_max = max.(xy_max, face_centroid[1:2])
                end
            end
        end
        Δ = (xy_max[1] - xy_min[1], xy_max[2] - xy_min[2], z_max - z_min)
        @assert all(x -> x > 0, Δ) "Cell dimensions were zero? Computed $Δ for cell $index."
    end
    return Δ
end

function compute_peaceman_index(Δ, K, radius;
        dir::Symbol = :z,
        net_to_gross = 1.0,
        constant = 0.14,
        Kh = nothing,
        drainage_radius = nothing,
        skin = 0,
        check = true
    )
    K_d = diag(K)
    Δx, Δy, Δz = Δ
    Δz *= net_to_gross
    if dir == :x || dir == :X
        L, d1, d2 = Δx, Δy, Δz
        i, j = 2, 3
    elseif dir == :y || dir == :Y
        d1, L, d2 = Δx, Δy, Δz
        i, j = 1, 3
    else
        d1, d2, L = Δx, Δy, Δz
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
    ke  = sqrt(k1*k2)
    if isnothing(drainage_radius)
        re1 = 2 * constant * sqrt((d1^2)*sqrt(k21) + (d2^2)*sqrt(k12))
        re2 = k21^(1/4) + k12^(1/4)
        re = kratio(re1, re2)
    else
        re = drainage_radius
    end

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

function Jutul.discretize_domain(d::DataDomain, system::Union{MultiPhaseSystem, CompositeSystem{:Reservoir, T}}, ::Val{:default}; kwarg...) where T
    return discretized_domain_tpfv_flow(d; kwarg...)
end


function discretized_domain_tpfv_flow(domain::Jutul.DataDomain;
        general_ad = false,
        kgrad = nothing,
        upwind = nothing
    )
    N = domain[:neighbors]
    g = physical_representation(domain)
    nc = number_of_cells(g)
    # defaulted_disc = isnothing(kgrad) && isnothing(upwind)
    kgrad_is_tpfa = (isnothing(kgrad) || eltype(kgrad) == TPFA || kgrad == :tpfa)
    upw_is_tpfa = (isnothing(upwind) || eltype(upwind) == SPU || upwind == :spu)

    is_tpfa = kgrad_is_tpfa && upw_is_tpfa
    if is_tpfa
        if general_ad
            d = PotentialFlow(N, nc)
        else
            d = TwoPointPotentialFlowHardCoded(N, nc)
        end
    else
        if kgrad == :tpfa_test
            # Fallback version - use generic FVM assembly with TPFA.
            kgrad = nothing
        end
        if kgrad isa Symbol
            if kgrad == :tpfa
                kgrad = nothing
            else
                K = domain[:permeability]
                g = UnstructuredMesh(g)
                T_base = reservoir_transmissibility(domain)
                kgrad = Jutul.NFVM.ntpfa_decompose_faces(g, K, kgrad, tpfa_trans = T_base)
            end
        else
            @assert isnothing(kgrad) || kgrad isa AbstractVector
        end
        if general_ad
            ad_flag = :generic
        else
            ad_flag = :fvm
        end
        d = PotentialFlow(N, nc, kgrad = kgrad, upwind = upwind, ad = ad_flag)
    end
    disc = (mass_flow = d, heat_flow = d)
    G = MinimalTPFATopology(N, ncells = nc)
    return DiscretizedDomain(G, disc)
end

function Jutul.discretize_domain(d::DataDomain{W}, system::Union{MultiPhaseSystem, CompositeSystem{:Reservoir, T}}, ::Val{:default}; kwarg...) where {W<:Union{SimpleWell, MultiSegmentWell}, T}
    return discretized_domain_well(physical_representation(d); kwarg...)
end

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
