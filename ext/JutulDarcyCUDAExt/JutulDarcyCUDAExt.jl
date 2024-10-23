module JutulDarcyCUDAExt
    using Jutul, JutulDarcy, CUDA, LinearAlgebra, SparseArrays
    include("ilu0.jl")
    include("krylov.jl")
    include("cuda_utils.jl")
end
