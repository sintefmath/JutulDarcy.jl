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
struct MinimalTPFATopology{N} <: ReservoirGrid
    nc::Int
    neighborship::N
end

function MinimalTPFATopology(N::T; ncells = maximum(N)) where T<:AbstractMatrix
    @assert size(N, 1) == 2
    @assert maximum(N) <= ncells
    @assert minimum(N) > 0
    return MinimalTPFATopology{T}(ncells, N)
end

Base.show(io::IO, g::MinimalTPFATopology) = print(io, "MinimalTPFATopology ($(number_of_cells(g)) cells, $(number_of_faces(g)) faces)")

function number_of_cells(G::MinimalTPFATopology)
    return G.nc
end

function declare_entities(G::MinimalTPFATopology)
    c = (entity = Cells(), count = number_of_cells(G)) # Cells equal to number of pore volumes
    f = (entity = Faces(), count = number_of_faces(G)) # Faces
    return [c, f]
end

function transfer(context::SingleCUDAContext, grid::MinimalTPFAGrid)
    pv = transfer(context, grid.pore_volumes)
    N = transfer(context, grid.neighborship)

    return MinimalTPFAGrid(pv, N)
end

function get_1d_reservoir(nc; L = 1.0, perm = 9.8692e-14, # 0.1 darcy
                         poro = 0.1, area = 1.0,
                         z_max = nothing)
    @assert nc > 1 "Must have at least two cells."
    g = CartesianMesh((nc, 1, 1), (L, sqrt(area), sqrt(area)))
    D = reservoir_domain(g, permeability = perm, porosity = poro)

    if !isnothing(z_max)
        # Manipulate z coordinates manually
        dz = z_max/nc
        z = (dz/2:dz:z_max-dz/2)'
        cc = D[:cell_centroids]
        for i in eachindex(z)
            cc[3, i] = z[i]
        end
    end
    return D
end

function subgrid_renumber_neighborship(g, cells, faces)
    N = g.neighborship
    nc = number_of_cells(g)
    renumeration = Dict{Int, Int}()
    for (i, c) in enumerate(cells)
        renumeration[c] = i
    end
    N_new = N[:, faces]
    for i in eachindex(N_new)
        old = N_new[i]
        @assert haskey(renumeration, old)
        new = renumeration[old]
        N_new[i] = new
    end
    return N_new
end

function Jutul.subgrid(g::MinimalTPFAGrid; cells = nothing, faces = nothing)
    pv = g.pore_volumes
    N_new = subgrid_renumber_neighborship(g, cells, faces)
    return MinimalTPFAGrid(pv[cells], N_new)
end

function Jutul.subgrid(g::MinimalTPFATopology; cells, faces)
    N_new = subgrid_renumber_neighborship(g, cells, faces)
    return MinimalTPFATopology(length(cells), N_new)
end
