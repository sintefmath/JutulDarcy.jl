function JutulDarcy.gpu_reduce_residual!(r_p, w_p, r)
    ncomp, np = size(w_p)
    @assert eltype(r_p) == eltype(w_p) == eltype(r)
    # @assert length(w_p) == length(r)
    @assert np == length(r_p)
    ncell = length(r) ÷ ncomp
    kernel = CUDA.@cuda launch = false reduce_residual_kernel!(r_p, w_p, r)
    threads, blocks = kernel_helper(kernel, r_p)
    CUDA.@sync begin
        kernel(r_p, w_p, r; threads, blocks)
    end
    if np > ncell
        @. r_p[ncell+1:end] = 0
    end
    return r_p
end

function reduce_residual_kernel!(r_p, w_p, r_ps)
    index = threadIdx().x
    stride = blockDim().x
    T = eltype(r_p)
    ncomp = size(w_p, 1)
    for cell in index:stride:(length(r_ps)÷ncomp)
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

function JutulDarcy.gpu_increment_pressure!(x, dp, bz)
    n = length(x)
    @assert bz > 0
    x_view = view(x, 1:bz:(n+1-bz))
    if length(dp)*bz > n
        dp = view(dp, 1:(n÷bz))
        @. x_view += dp
    else
        @. x_view += dp
    end
    return x
end

function kernel_helper(kernel, V)
    config = launch_configuration(kernel.fun)
    threads = min(length(V), config.threads)
    blocks = cld(length(V), threads)

    return (threads, blocks)
end