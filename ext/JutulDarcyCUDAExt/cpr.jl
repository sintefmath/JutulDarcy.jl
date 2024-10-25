function JutulDarcy.gpu_reduce_residual!(r_p, w_p, r)
    # TODO: Actually use CUDA threads here.
    @assert eltype(r_p) == eltype(w_p) == eltype(r)
    CUDA.@cuda reduce_residual_kernel!(r_p, w_p, r)
    return r_p
end

function reduce_residual_kernel!(r_p, w_p, r_full)
    index = threadIdx().x
    stride = blockDim().x
    T = eltype(r_p)
    ncomp = size(w_p, 1)
    for cell in index:stride:length(r_p)
        v = zero(T)
        for comp in axes(w_p, 1)
            v += w_p[comp, cell] * r_full[(cell-1)*ncomp + comp]
        end
        r_p[cell] = v
    end
    return nothing
end

function JutulDarcy.gpu_increment_pressure!(x, dp)
    n = length(x)
    bz = div(n, length(dp))
    @assert bz > 0
    CUDA.@cuda reduce_residual_kernel!(x, dp, bz)
    return x
end

function gpu_increment_pressure_kernel!(x, dp, bz)
    index = threadIdx().x
    stride = blockDim().x
    for cell in index:stride:length(dp)
        x[(cell-1)*bz + 1] = dp[cell]
    end
    return nothing
end
