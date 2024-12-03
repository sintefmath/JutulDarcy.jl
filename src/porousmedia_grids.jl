
import Base.eltype

function number_of_cells(G::ReservoirGrid)
    return length(G.pore_volumes)
end

struct MinimalTPFATopology{N} <: ReservoirGrid
    nc::Int
    neighborship::N
end

function MinimalTPFATopology(N::T; ncells = maximum(N, init = 1)) where T<:AbstractMatrix
    @assert size(N, 1) == 2
    @assert maximum(N, init = 1) <= ncells
    @assert minimum(N, init = 1) > 0
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

"""
    get_1d_reservoir(nc;
        L = 1.0,
        perm = 9.8692e-14, # 0.1 darcy
        poro = 0.1,
        area = 1.0,
        z_max = nothing
    )

Utility function for setting up a 1D reservoir domain with `nc` cells and length
`L`. The [`reservoir_domain`](@ref) function is generally preferred and this
function is kept for backwards compatibility.
"""
function get_1d_reservoir(nc;
        L = 1.0,
        perm = 9.8692e-14, # 0.1 darcy
        poro = 0.1,
        area = 1.0,
        z_max = nothing
    )
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

function Jutul.subgrid(g::MinimalTPFATopology; cells, faces)
    N_new = subgrid_renumber_neighborship(g, cells, faces)
    return MinimalTPFATopology(length(cells), N_new)
end
