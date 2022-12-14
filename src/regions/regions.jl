function check_regions(regions::AbstractVector)
    @assert minimum(regions) > 0
    @assert eltype(regions)<:Integer
end

function check_regions(regions::Nothing)
    # Ok.
end

@inline function region(pv::DeckPhaseVariables, cell)
    return region(pv.regions, cell)
end

@inline function region(r::AbstractVector, cell)
    return @inbounds r[cell]
end

@inline function region(::Nothing, cell)
    return 1
end
@inline number_of_regions(regions::Nothing) = 1
@inline number_of_regions(regions::AbstractVector) = 1

@inline function table_by_region(tab, reg)
    return tab[reg]
end

function region_wrap(x::Nothing, regions)
    return x
end

function region_wrap(x::Tuple, regions::Nothing)
    @assert length(x) == 1
    return x
end

function region_wrap(x::Tuple, regions::AbstractArray)
    @assert length(x) >= maximum(regions)
    return x
end

function region_wrap(x::AbstractVector, regions)
    return region_wrap(Tuple(x), regions)
end

function region_wrap(x::Any, regions)
    return region_wrap((x,), regions)
end
