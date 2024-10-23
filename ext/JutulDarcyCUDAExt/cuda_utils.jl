function JutulDarcy.build_gpu_block_system(Ti, Tv, sz::Tuple{Int, Int}, blockDim::Int, rowptr, colval, nzval, r0)
    rowPtr = CuVector{Ti}(rowptr)
    colVal = CuVector{Ti}(colval)
    nzVal = CuVector{Tv}(nzval)
    dims = blockDim.*sz
    dir = 'C'
    nnz = length(nzVal)÷(blockDim^2)
    J_bsr = CUDA.CUSPARSE.CuSparseMatrixBSR{Tv, Ti}(rowPtr, colVal, nzVal, dims, blockDim, dir, nnz)
    r_cu = CuVector{Tv}(vec(r0))
    return (J_bsr, r_cu)
end

function JutulDarcy.update_gpu_block_system!(J, blockDim, nzval, ϵ = 1e-12)
    if ϵ > 0.0
        for (i, v) in enumerate(nzval)
            if abs(v) < ϵ
                if v > 0
                    v = ϵ
                else
                    v = ϵ
                end
                nzval[i] = v
            end
        end
    end
    copyto!(J.nzVal, nzval)
end

function JutulDarcy.update_gpu_block_residual!(r_cu, blockDim, r)
    copyto!(r_cu, r)
end
