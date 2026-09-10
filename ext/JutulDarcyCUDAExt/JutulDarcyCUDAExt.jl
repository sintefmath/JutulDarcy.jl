module JutulDarcyCUDAExt
    using Jutul, JutulDarcy, CUDA, LinearAlgebra, SparseArrays
    import Jutul: @tic

    JutulDarcy.kernel_abstractions_backend(::Val{:ka_cuda}) = CUDA.CUDABackend()

    timeit_debug_enabled() = Jutul.timeit_debug_enabled()

    include("ilu0.jl")
    include("krylov.jl")
    include("cuda_utils.jl")
    include("cpr.jl")
end
