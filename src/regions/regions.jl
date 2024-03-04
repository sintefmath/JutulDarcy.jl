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

Base.@propagate_inbounds @inline function evaluate_table_by_region(tab, reg, arg...)
    return tab[reg](arg...)
end

Base.@propagate_inbounds @inline function evaluate_table_by_region(tab, ::Nothing, arg...)
    return tab(arg...)
end

Base.@propagate_inbounds @inline function table_by_region(tab, reg)
    return tab[reg]
end

@inline function table_by_region(tab::Nothing, reg)
    return tab
end

function region_wrap(x::Nothing, regions = missing)
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

function region_wrap(x::Tuple, regions::Missing)
    return x
end

function region_wrap(x::AbstractVector, regions = missing)
    return region_wrap(Tuple(x), regions)
end

"""
    region_wrap(x::Any, regions = missing)

Turn one or more functions/tables into a Tuple for access using
`table_by_region`. Single tables will be turned into a single-element tuple.
"""
function region_wrap(x::Any, regions = missing)
    return region_wrap((x,), regions)
end
