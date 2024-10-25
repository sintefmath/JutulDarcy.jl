module JutulDarcyAMGXExt
    using Jutul, JutulDarcy, CUDA, LinearAlgebra, SparseArrays, AMGX
    import Jutul: @tic

    timeit_debug_enabled() = Jutul.timeit_debug_enabled()

    function JutulDarcy.gpu_update_preconditioner!(prec::CPRPreconditioner{JutulDarcy.AMGXPreconditioner, <:Any}, sys, model, storage, recorder, executor, krylov, J_bsr, r_cu)
        error()
        Jutul.update_preconditioner!(krylov.preconditioner, J_bsr, r_cu, model.context, executor)
    end

end
