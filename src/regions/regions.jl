function check_regions(regions)
    if !isnothing(regions)
        regions<:AbstractVector
        @assert minimum(regions) > 0
    end
end


@inline region(pv::DeckPhaseVariables, cell) = region(pv.regions, cell)
@inline region(r::AbstractVector, cell) = @inbounds r[cell]
@inline region(::Nothing, cell) = 1

@inline number_of_regions(regions::Nothing) = 1
@inline number_of_regions(regions::AbstractVector) = 1

@inline function table_by_region(tab, reg)
    return tab[reg]
end

function region_wrap(x::Tuple, regions::Nothing)
    @assert length(x) == 1
    return x
end

function region_wrap(x::Tuple, regions::AbstractArray)
    @assert length(x) <= maximum(regions)
    return x
end

function region_wrap(x::AbstractVector, regions)
    return region_wrap(Tuple(x...), regions)
end

function region_wrap(x::Any, regions)
    return region_wrap((x,), regions)
end
