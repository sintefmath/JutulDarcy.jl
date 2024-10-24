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

function JutulDarcy.build_gpu_schur_system(Ti, Tv, bz, lsys::LinearizedSystem)
    return nothing
end

function JutulDarcy.build_gpu_schur_system(Ti, Tv, bz, lsys::MultiLinearizedSystem)
    @assert size(lsys.subsystems) == (2, 2)
    @assert length(lsys.schur_buffer) == 2
    buf_full_cpu = lsys.schur_buffer[1]
    buf_1_cpu, buf_2_cpu = lsys.schur_buffer[2]
    # B_cpu = lsys[1, 1].jac # Already transferred, not needed
    C_cpu = lsys[1, 2].jac
    D_cpu = lsys[2, 1].jac
    # E_cpu = lsys[2, 2].jac # Only factorization needed I think.

    E_factor = only(lsys.factor.factor)
    E_L_cpu = E_factor.L
    E_U_cpu = E_factor.U

    to_gpu(x) = copy_to_gpu(x, Tv, Ti)

    D = to_gpu(D_cpu)
    E_L = to_gpu(E_L_cpu)
    E_U = to_gpu(E_U_cpu)
    C = to_gpu(C_cpu)

    buf_1 = to_gpu(buf_1_cpu)
    buf_2 = to_gpu(buf_2_cpu)
    # Operations and format:
    # B C
    # D E
    # mul!(b_buf_2, D_i, x)
    # ldiv!(b_buf_1, E_i, b_buf_2)
    # mul!(res, C_i, b_buf_1, -α, true)

    return Dict(
        :C => C,
        :D => D,
        :buf_1 => buf_1,
        :buf_2 => buf_2,
        :E_L => E_L,
        :E_U => E_U
    )
end

function JutulDarcy.update_gpu_schur_system!(schur::Nothing, lsys)
    return schur
end

function JutulDarcy.update_gpu_schur_system!(schur, lsys)
    C_cpu = lsys[1, 2].jac
    copyto!(schur[:C].nzVal, C_cpu.nzval)
    D_cpu = lsys[2, 1].jac
    copyto!(schur[:D].nzVal, D_cpu.nzval)
    E_factor = only(lsys.factor.factor)
    copyto!(schur[:E_L].nzVal, E_factor.L.nzval)
    copyto!(schur[:E_U].nzVal, E_factor.U.nzval)
    return schur
end

function copy_to_gpu(x, Ti, Tv)
    error("Internal converter not set up for $(typeof(x))")
end

function copy_to_gpu(x::SparseMatrixCSC{Tvc, Tic}, Tv, Ti) where {Tvc, Tic}
    if Tvc != Tv || Tic != Ti
        x = SparseMatrixCSC{Tv, Ti}(x)
    end
    return CUDA.CUSPARSE.CuSparseMatrixCSC(x)
end

function copy_to_gpu(x::Vector{Tvc}, Tv, Ti) where {Tvc}
    return convert(CuVector{Tv}, x)
end
