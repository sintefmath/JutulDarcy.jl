export MinimalTPFAGrid
export subgrid

export transfer, get_1d_reservoir

import Base.eltype

# TPFA grid
"""
    MinimalTPFAGrid(ϕ, N)

Generate a minimal grid suitable only for two-point flux discretization (TPFA) for given
pore-volumes `ϕ` and a neighborship matrix `N` with size `(2, n)` where `n` is the number
of internal faces.
"""
struct MinimalTPFAGrid{V, N} <: ReservoirGrid
    pore_volumes::V
    trans::V
    gdz::V
    neighborship::N
    function MinimalTPFAGrid(pv::V, N::M; trans::V = ones(size(N, 2)), gdz::V = zeros(size(N, 2))) where {V<:AbstractVector, M<:AbstractMatrix}
        nc = length(pv)
        pv::AbstractVector
        @assert size(N, 1) == 2  "Two neighbors per face"
        if length(N) > 0
            @assert minimum(N) > 0   "Neighborship entries must be positive."
            @assert maximum(N) <= nc "Neighborship must be limited to number of cells."
        end
        @assert all(pv .> 0)     "Pore volumes must be positive"
        new{V, M}(pv, trans, gdz, N)
    end
end
Base.show(io::IO, g::MinimalTPFAGrid) = print(io, "MinimalTPFAGrid ($(number_of_cells(g)) cells, $(number_of_faces(g)) faces)")


function number_of_cells(G::ReservoirGrid)
    return length(G.pore_volumes)
end

function declare_entities(G::MinimalTPFAGrid)
    c = (entity = Cells(), count = number_of_cells(G)) # Cells equal to number of pore volumes
    f = (entity = Faces(), count = number_of_faces(G)) # Faces
    return [c, f]
end

function transfer(context::SingleCUDAContext, grid::MinimalTPFAGrid)
    pv = transfer(context, grid.pore_volumes)
    N = transfer(context, grid.neighborship)

    return MinimalTPFAGrid(pv, N)
end

function get_1d_reservoir(nc; L = 1, perm = 9.8692e-14, # 0.1 darcy
                         poro = 0.1, area = 1,
                         z_max = nothing, general_ad = false)
    @assert nc > 1 "Must have at least two cells."
    nf = nc-1
    N = vcat((1:nf)', (2:nc)')
    dx = L/nc
    cell_centroids = vcat((dx/2:dx:L-dx/2)', ones(2, nc))
    face_centroids = vcat((dx:dx:L-dx)', ones(2, nf))
    face_areas = ones(nf)
    face_normals = vcat(ones(1, nf), zeros(2, nf))

    function expand(x, nc)
        repeat([x], nc)
    end
    function expand(x::AbstractVector, nc)
        x
    end
    perm = expand(perm, nc)
    if isa(perm, AbstractVector)
        perm = copy(perm')
    end
    volumes = repeat([area.*dx], nc)

    pv = poro.*volumes
    nc = length(pv)

    @debug "Data unpack complete. Starting transmissibility calculations."
    # Deal with face data
    if isnothing(z_max)
        z = zeros(nc)
    else
        dz = z_max/nc
        z = (dz/2:dz:z_max-dz/2)'
    end
    g = gravity_constant
    T_hf = compute_half_face_trans(cell_centroids, face_centroids, face_normals, face_areas, perm, N)
    T = compute_face_trans(T_hf, N)
    gdz = compute_face_gdz(N, z)
    G = MinimalTPFAGrid(pv, N, trans = T, gdz = gdz)
    if general_ad
        flow = PotentialFlow(N, nc)
    else
        flow = TwoPointPotentialFlowHardCoded(G, ncells = nc)
    end
    disc = (mass_flow = flow, heat_flow = flow)
    D = DiscretizedDomain(G, disc)
    return D
end

function Jutul.subgrid(g::MinimalTPFAGrid; cells = nothing, faces = nothing)
    pv, N = g.pore_volumes, g.neighborship
    nc = number_of_cells(g)

    renumeration = zeros(Integer, nc)
    for (i, c) in enumerate(cells)
        renumeration[c] = i
    end
    N_new = N[:, faces]
    for i in eachindex(N_new)
        old = N_new[i]
        new = renumeration[old]
        @assert new != 0
        N_new[i] = new
    end
    return MinimalTPFAGrid(pv[cells], N_new)
end
