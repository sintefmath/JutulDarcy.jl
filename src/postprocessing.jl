# Generic fallback
delta_state(xa, xb, relative=false) = missing

# Scalars
function delta_state(xa::Number, xb::Number, relative=false)
    relative ? (xa - xb) / xb : (xa - xb)
end

# Arrays (matrix/tensor/etc.; vectors handled by more specific methods below)
function delta_state(xa::AbstractArray, xb::AbstractArray, relative=false)
    size(xa) == size(xb) || return missing
    if eltype(xa) <: Number && eltype(xb) <: Number
        println("Should have been here!")
        return relative ? (xa .- xb) ./ xb : (xa .- xb)
    else
        println("Should not have been here!")
        return map((a, b) -> delta_state(a, b, relative), xa, xb)
    end
end

# Dict/OrderedDict (and other AbstractDict)
function delta_state(xa::AbstractDict, xb::AbstractDict, relative=false)
    # dx = OrderedDict{Any, Any}()
    dx = deepcopy(xa)
    for (k, va) in xa
        println("Processing key: $k")
        dx[k] = haskey(xb, k) ? delta_state(va, xb[k], relative) : missing
        println("Type of dx[$k]: $(typeof(dx[k]))")
    end
    return dx
end

# Vector vs vector (e.g., vector of states)
function delta_state(xa::AbstractVector, xb::AbstractVector, relative=false)
    length(xa) == length(xb) || return missing
    dx = Vector{Any}(undef, length(xa))
    @inbounds for i in eachindex(xa, xb)
        dx[i] = delta_state(xa[i], xb[i], relative)
    end
    return dx
end

# Vector of states vs single reference state (dict/ordered dict)
function delta_state(xa::AbstractVector, xb::AbstractDict, relative=false)
    dx = Vector{Any}(undef, length(xa))
    @inbounds for i in eachindex(xa)
        dx[i] = delta_state(xa[i], xb, relative)
    end
    return dx
end

# Single state vs vector of reference states
function delta_state(xa::AbstractDict, xb::AbstractVector, relative=false)
    dx = Vector{Any}(undef, length(xb))
    @inbounds for i in eachindex(xb)
        dx[i] = delta_state(xa, xb[i], relative)
    end
    return dx
end