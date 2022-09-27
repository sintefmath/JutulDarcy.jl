
expand_to_phases(v::Real, nph) = SVector{nph}([v for i in 1:nph])
function expand_to_phases(v::AbstractVector, nph)
    @assert length(v) == nph
    return SVector{nph}(v)
end
