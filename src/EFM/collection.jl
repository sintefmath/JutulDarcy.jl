"""
    FractureCollection(fractures)

Group a fracture vector by concrete type for specialized preprocessing. Input
vectors may have any element type, including `Any`, and may be views. Every
entry must be an `AbstractEmbeddedFracture`. Each vector in `groups` has a
concrete element type, and the corresponding vector in `indices` contains its
original fracture IDs. Indexing the collection retains the input order.

The collection copies the input vector's entries into its groups. Treat the
group vectors as read-only; construct a new collection after adding fractures.
Both EFM and DFM preprocessing create a collection automatically.
"""
struct FractureCollection{G <: Tuple, I <: Tuple} <: AbstractVector{AbstractEmbeddedFracture}
    groups::G
    indices::I
    lookup::Vector{Tuple{Int, Int}}
end

FractureCollection(fractures::FractureCollection) = fractures
FractureCollection(fracture::AbstractEmbeddedFracture) = FractureCollection((fracture,))

function FractureCollection(fractures::Union{AbstractVector, Tuple})
    groups = Any[]
    indices = Vector{Int}[]
    lookup = Tuple{Int, Int}[]
    types = Dict{DataType, Int}()
    for (id, fracture) in enumerate(fractures)
        fracture isa AbstractEmbeddedFracture || throw(ArgumentError("fracture at position $id must be an AbstractEmbeddedFracture, got $(typeof(fracture))"))
        group = get(types, typeof(fracture), 0)
        if group == 0
            push!(groups, Vector{typeof(fracture)}())
            push!(indices, Int[])
            group = length(groups)
            types[typeof(fracture)] = group
        end
        push!(groups[group], fracture)
        push!(indices[group], id)
        push!(lookup, (group, length(groups[group])))
    end
    return FractureCollection(Tuple(groups), Tuple(indices), lookup)
end

Base.size(fractures::FractureCollection) = (length(fractures.lookup),)
Base.IndexStyle(::Type{<:FractureCollection}) = IndexLinear()

function Base.getindex(fractures::FractureCollection, id::Int)
    @boundscheck checkbounds(fractures, id)
    group, local_index = fractures.lookup[id]
    return fractures.groups[group][local_index]
end

function fracture_property_defaults(fractures::FractureCollection, ::Type{T}) where T
    n = length(fractures)
    values = (aperture = Vector{T}(undef, n), permeability = Vector{T}(undef, n),
        porosity = Vector{T}(undef, n))
    foreach(fractures.groups, fractures.indices) do group, ids
        fracture_property_defaults!(values, group, ids)
    end
    return values
end

function fracture_property_defaults!(values, group, ids)
    for (fracture, id) in zip(group, ids)
        properties = fracture.properties
        values.aperture[id] = properties.aperture
        values.permeability[id] = properties.permeability
        values.porosity[id] = properties.porosity
    end
    return values
end
