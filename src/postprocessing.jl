function delta_state(x::AbstractArray, x0::AbstractArray, relative=false)
    if length(x) != length(x0)
        return missing
    else
        if relative
            return (x .- x0) ./ x0
        else
            return x .- x0
        end
    end
end

function delta_state(x, x0::Dict, relative=false)
    dx = deepcopy(x)
    for (k, v) in x
        if haskey(x0, k)
            dx[k] = delta_state(x[k], x0[k], relative)
        else
            dx[k] = missing
        end
    end
    return dx
end

function delta_state(x::Vector, x0::Dict, relative=false)
    dx = deepcopy(x)
    for (n, xn) in enumerate(x)
        dx[n] = delta_state(xn, x0, relative)
    end
    return dx
end