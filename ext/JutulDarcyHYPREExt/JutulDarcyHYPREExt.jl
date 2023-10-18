module JutulDarcyHYPREExt
    using JutulDarcy
    using HYPRE
    using Jutul
    using SparseArrays
    using PrecompileTools
    include("cpr.jl")

    # @compile_workload begin
    #     targets = [(true, :csc), (true, :csr)]
    #     # MPI, trivial partition
    #     JutulDarcy.precompile_darcy_multimodels(targets,
    #         dims = (4, 1, 1),
    #         precond = :cpr,
    #         amg_type = :hypre
    #     )
    # end
end
