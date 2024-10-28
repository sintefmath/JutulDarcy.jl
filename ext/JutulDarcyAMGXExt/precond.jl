function JutulDarcy.gpu_update_preconditioner!(prec::JutulDarcy.AMGXPreconditioner, lsys, model, storage, recorder, executor, krylov, J_bsr, r_cu, op)
    nzval = nonzeros(J_bsr)
    bz = Jutul.block_size(lsys[1, 1])
    if haskey(prec.data, :storage)
        s = prec.data[:storage]
        AMGX.replace_coefficients!(s.matrix, nzval)
    else
        Tv = eltype(J_bsr)
        if Tv == Float64
            amgx_mode = AMGX.dDDI
        else
            amgx_mode = AMGX.dFFI
        end
        n, m = size(J_bsr)
        @assert n == m
        N = n ÷ bz
        config = AMGX.Config(prec.settings)
        s = AMGXStorage(config, amgx_mode)

        AMGX.set_zero!(s.x, N, block_dim = bz)
        AMGX.set_zero!(s.r, N, block_dim = bz)

        row_ptr = J_bsr.rowPtr
        colval = J_bsr.colVal
        AMGX.upload!(
            s.matrix,
            row_ptr,
            colval,
            nzval,
            block_dims = (bz, bz)
        )
        prec.data[:storage] = s
        prec.data[:n] = N
        prec.data[:block_size] = bz
    end
    AMGX.setup!(s.solver, s.matrix)
end

function Jutul.operator_nrows(prec::JutulDarcy.AMGXPreconditioner)
    return prec.data[:n]*prec.data[:block_size]
end

function Jutul.apply!(x, amgx::JutulDarcy.AMGXPreconditioner, r)
    # TODO: This sync call is maybe needd?
    AMGX.CUDA.synchronize()
    s = amgx.data[:storage]
    bz = amgx.data[:block_size]::Int
    n = amgx.data[:n]::Int
    AMGX.upload!(s.r, r, block_dim = bz)
    AMGX.set_zero!(s.x, n, block_dim = bz)
    AMGX.solve!(s.x, s.solver, s.r)
    AMGX.download!(x, s.x)
end
