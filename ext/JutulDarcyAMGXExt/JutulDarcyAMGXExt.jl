module JutulDarcyAMGXExt
    using JutulDarcy, Jutul, AMGX, AMGX.CUDA, LinearAlgebra, SparseArrays
    import JutulDarcy: AMGXPreconditioner
    import Jutul: @tic

    timeit_debug_enabled() = Jutul.timeit_debug_enabled()

    include("cpr.jl")

    function __init__()
        # TODO: Figure out a way to not always do this?
        AMGX.initialize()
    end
end
