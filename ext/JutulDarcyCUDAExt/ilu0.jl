function ilu0_invert!(y, A::CUDA.CUSPARSE.CuSparseMatrixBSR, Tv, c_char, data)
    if Tv == Float32
        bname = CUDA.CUSPARSE.cusparseSbsrsv2_bufferSize
        aname = CUDA.CUSPARSE.cusparseSbsrsv2_analysis
        sname = CUDA.CUSPARSE.cusparseSbsrsv2_solve
    else
        bname = CUDA.CUSPARSE.cusparseDbsrsv2_bufferSize
        aname = CUDA.CUSPARSE.cusparseDbsrsv2_analysis
        sname = CUDA.CUSPARSE.cusparseDbsrsv2_solve
    end
    if c_char == 'U'
        d_char = 'N'
        kval = :ilu_upper
    else
        @assert c_char == 'L'
        d_char = 'U'
        kval = :ilu_lower
    end
    transa = 'N'
    uplo = c_char
    diag = d_char
    alpha = 1.0
    index = 'O'
    X = y
    m,n = size(A)
    if m != n
        throw(DimensionMismatch("A must be square, but has dimensions ($m,$n)!"))
    end

    mb = cld(m, A.blockDim)

    is_analysed = haskey(data, kval)
    # is_analysed = false
    if is_analysed
        info, desc = data[kval]
    else
        info = CUDA.CUSPARSE.bsrsv2Info_t[0]
        desc = CUDA.CUSPARSE.CuMatrixDescriptor('G', uplo, diag, index)
        CUDA.CUSPARSE.cusparseCreateBsrsv2Info(info)
        data[kval] = (info, desc)
    end

    function bufferSize()
        out = Ref{Cint}(1)
        bname(CUDA.CUSPARSE.handle(), A.dir, transa, mb, A.nnzb,
            desc, nonzeros(A), A.rowPtr, A.colVal, A.blockDim,
            info[1], out)
        return out[]
    end
    CUDA.CUSPARSE.with_workspace(bufferSize) do buffer
        if !is_analysed
            aname(CUDA.CUSPARSE.handle(), A.dir, transa, mb, A.nnzb,
                    desc, nonzeros(A), A.rowPtr, A.colVal, A.blockDim,
                    info[1], CUDA.CUSPARSE.CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer)
            posit = Ref{Cint}(1)
            CUDA.CUSPARSE.cusparseXbsrsv2_zeroPivot(CUDA.CUSPARSE.handle(), info[1], posit)
            if posit[] >= 0
                error("Structural/numerical zero in A at ($(posit[]),$(posit[])))")
            end
        end
        sname(CUDA.CUSPARSE.handle(), A.dir, transa, mb, A.nnzb,
                alpha, desc, nonzeros(A), A.rowPtr, A.colVal,
                A.blockDim, info[1], X, X,
                CUDA.CUSPARSE.CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer)
    end
    # Alternatively we could not cache setup and do the following:
    # if c_char == 'U'
    #     CUDA.CUSPARSE.sv2!('N', 'U', 'N', 1.0, A, y, 'O')
    # else
    #     @assert c_char == 'L'
    #     CUDA.CUSPARSE.sv2!('N', 'L', 'U', 1.0, A, y, 'O')
    # end
    return y
end

function ilu0_gpu_update!(A, Tv, data)
    if Tv == Float32
        bname = CUDA.CUSPARSE.cusparseSbsrilu02_bufferSize
        aname = CUDA.CUSPARSE.cusparseSbsrilu02_analysis
        sname = CUDA.CUSPARSE.cusparseSbsrilu02
    else
        bname = CUDA.CUSPARSE.cusparseDbsrilu02_bufferSize
        aname = CUDA.CUSPARSE.cusparseDbsrilu02_analysis
        sname = CUDA.CUSPARSE.cusparseDbsrilu02
    end
    mb = div(size(A, 1), A.blockDim)
    ilu0_is_analyzed = haskey(data, :ilu0)
    if ilu0_is_analyzed
        info, desc = data[:ilu0]
    else
        info = CUDA.CUSPARSE.ILU0InfoBSR()
        desc = CUDA.CUSPARSE.CuMatrixDescriptor('G', 'U', 'N', 'O')
        data[:ilu0] = (info, desc)
    end

    function bufferSize()
        out = Ref{Cint}(1)
        bname(CUDA.CUSPARSE.handle(), A.dir, mb, nnz(A), desc, nonzeros(A),
               A.rowPtr, A.colVal, A.blockDim, info, out)
        return out[]
    end

    policy = CUDA.CUSPARSE.CUSPARSE_SOLVE_POLICY_USE_LEVEL
    CUDA.CUSPARSE.with_workspace(bufferSize) do buffer
        if !ilu0_is_analyzed
            aname(
                CUDA.CUSPARSE.handle(),
                A.dir,
                mb,
                nnz(A),
                desc,
                nonzeros(A),
                A.rowPtr,
                A.colVal,
                A.blockDim,
                info,
                policy,
                buffer
            )
            posit = Ref{Cint}(1)
            CUDA.CUSPARSE.cusparseXbsrilu02_zeroPivot(CUDA.CUSPARSE.handle(), info, posit)
            if posit[] >= 0
                error("Structural/numerical zero in A at ($(posit[]),$(posit[])))")
            end
        end
        sname(
            CUDA.CUSPARSE.handle(),
            A.dir,
            mb,
            SparseArrays.nnz(A),
            desc,
            nonzeros(A),
            A.rowPtr,
            A.colVal,
            A.blockDim,
            info,
            policy,
            buffer
        )
    end
end


function Jutul.update_preconditioner!(ilu::ILUZeroPreconditioner, J_bsr::CUDA.CUSPARSE.CuSparseMatrixBSR{Tv, Ti}, b, context, executor) where {Tv, Ti}
    if isnothing(ilu.factor)
        data = Dict()
        data[:ilu0_gpu] = copy(J_bsr)
        ilu.dim = size(J_bsr)
        ilu.factor = data
    end
    data = ilu.factor
    ilu_gpu = data[:ilu0_gpu]
    copyto!(ilu_gpu.nzVal, J_bsr.nzVal)
    ilu0_gpu_update!(ilu_gpu, Tv, data)
end

function Jutul.apply!(y::CuVector{Tv}, f::ILUZeroPreconditioner, x::CuVector{Tv}) where Tv
    copyto!(y, x)
    # Equivialent but does not skip analyse
    # sv2!('N', 'L', 'U', 1.0, P, y, 'O')
    # sv2!('N', 'U', 'N', 1.0, P, y, 'O')
    P = f.factor[:ilu0_gpu]
    ilu0_invert!(y, P, Tv, 'L', f.factor)
    ilu0_invert!(y, P, Tv, 'U', f.factor)
    CUDA.synchronize()
end
