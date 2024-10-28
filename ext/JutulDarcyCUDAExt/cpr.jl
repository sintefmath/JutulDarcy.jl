function JutulDarcy.gpu_reduce_residual!(r_p, w_p, r)
    # TODO: Actually use CUDA threads here.
    @assert eltype(r_p) == eltype(w_p) == eltype(r)
    @assert length(w_p) == length(r)
    @assert size(w_p, 2) == length(r_p)
    CUDA.@cuda reduce_residual_kernel!(r_p, w_p, r)
    return r_p
end

function reduce_residual_kernel!(r_p, w_p, r_ps)
    index = threadIdx().x
    stride = blockDim().x
    T = eltype(r_p)
    ncomp = size(w_p, 1)
    for cell in index:stride:length(r_p)
        v = zero(T)
        for comp in 1:ncomp
            ix = (cell-1)*ncomp + comp
            @inbounds wi = w_p[comp, cell]
            @inbounds ri = r_ps[ix]
            v += wi*ri
        end
        @inbounds r_p[cell] = v
    end
    return nothing
end

function JutulDarcy.gpu_increment_pressure!(x, dp)
    n = length(x)
    bz = div(n, length(dp))
    @assert bz > 0
    CUDA.@sync CUDA.@cuda gpu_increment_pressure_kernel!(x, dp, bz)
    # kernel = CUDA.@cuda launch=false gpu_increment_pressure_kernel!(x, dp, bz)
    # threads, blocks = kernel_helper(kernel, x)
    # CUDA.@sync begin
    #    kernel(x, dp, bz; threads, blocks)
    # end
    return x
end

function gpu_increment_pressure_kernel!(x, dp, bz)
    index = threadIdx().x
    stride = blockDim().x
    for cell in index:stride:length(dp)
        @inbounds x[(cell-1)*bz + 1] += dp[cell]
    end
    return nothing
end

function kernel_helper(kernel, V)
    config = launch_configuration(kernel.fun)
    threads = min(length(V), config.threads)
    blocks = cld(length(V), threads)

    return (threads, blocks)
end