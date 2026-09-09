module JutulDarcyCUDAExt
    using Jutul, JutulDarcy, CUDA, LinearAlgebra, SparseArrays
    import Jutul: @tic

    function __init__()
        JutulDarcy.register_kernel_abstractions_backend!(
            :ka_cuda, () -> CUDA.CUDABackend())
    end

    timeit_debug_enabled() = Jutul.timeit_debug_enabled()

    include("ilu0.jl")
    include("krylov.jl")
    include("cuda_utils.jl")
    include("cpr.jl")
end
