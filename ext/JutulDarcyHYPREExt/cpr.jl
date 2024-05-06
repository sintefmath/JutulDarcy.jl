function JutulDarcy.update_pressure_system!(A_p::HYPRE.HYPREMatrix, p_prec, A, w_p, ctx, executor)
    D = p_prec.data
    if !haskey(D, :assembly_helper)
        (; ilower, iupper) = A_p
        D[:assembly_helper] = Jutul.generate_hypre_assembly_helper(A, executor, ilower, iupper)
    end
    helper = D[:assembly_helper]
    I_buf, J_buf, V_buffers, = D[:assembly_helper]
    is_adjoint = Val(Jutul.represented_as_adjoint(matrix_layout(ctx)))
    ncomp = Val(size(w_p, 1))
    update_pressure_system_hypre!(I_buf, J_buf, V_buffers, A_p, A, w_p, executor, helper.n, is_adjoint, ncomp)
end

function update_pressure_system_hypre!(single_buf, longer_buf, V_buffers, A_p, A, w_p, executor, n, is_adjoint, ncomp)
    @assert length(single_buf) == 1
    (; iupper, ilower) = A_p
    @assert n == iupper - ilower + 1 "$(n-1) != $ilower -> $iupper"
    assembler = HYPRE.start_assemble!(A_p)
    assemble_into_hypre_psystem!(A_p, A, assembler, w_p, single_buf, longer_buf, V_buffers, executor, n, is_adjoint, ncomp)
    HYPRE.finish_assemble!(assembler)
end

function assemble_into_hypre_psystem!(A_p::HYPRE.HYPREMatrix, A::Jutul.StaticSparsityMatrixCSR, assembler, w_p, single_buf, longer_buf, V_buffers, executor, n, ::Val{is_adjoint}, ::Val{ncomp}) where {is_adjoint, ncomp}
    nzval = SparseArrays.nonzeros(A)
    cols = Jutul.colvals(A)
    for row in 1:n
        pos_ix = nzrange(A, row)
        k = length(pos_ix)
        I = single_buf
        I[1] = Jutul.executor_index_to_global(executor, row, :row)
        J = longer_buf
        resize!(J, k)
        V_buf = V_buffers[k]
        num_added = 0
        @inbounds for ki in 1:k
            ri = pos_ix[ki]
            col = cols[ri]
            A_block = nzval[ri]
            V_buf[ki] = JutulDarcy.reduce_to_pressure(A_block, w_p, row, ncomp, is_adjoint)
            J[ki] = Jutul.executor_index_to_global(executor, col, :column)
            num_added += 1
        end
        HYPRE.assemble!(assembler, I, J, V_buf)
    end
end

function assemble_into_hypre_psystem!(A_p::HYPRE.HYPREMatrix, A::SparseMatrixCSC, assembler, w_p, single_buf, longer_buf, V_buffers, executor, n, ::Val{is_adjoint}, ::Val{ncomp}) where {is_adjoint, ncomp}
    nzval = SparseArrays.nonzeros(A)
    rows = rowvals(A)
    for col in 1:n
        pos_ix = nzrange(A, col)
        k = length(pos_ix)
        J = single_buf
        J[1] = Jutul.executor_index_to_global(executor, col, :column)
        I = longer_buf
        resize!(I, k)
        V_buf = V_buffers[k]
        num_added = 0
        @inbounds for ki in 1:k
            ri = pos_ix[ki]
            row = rows[ri]
            A_block = nzval[ri]
            V_buf[ki] = JutulDarcy.reduce_to_pressure(A_block, w_p, row, ncomp, is_adjoint)
            I[ki] = Jutul.executor_index_to_global(executor, row, :row)
            num_added += 1
        end
        HYPRE.assemble!(assembler, I, J, V_buf)
    end
end

function JutulDarcy.create_pressure_system(p_prec::BoomerAMGPreconditioner, J, n)
    tmp = HYPRE.HYPREVector(zeros(3))
    comm = HYPRE.Internals.get_comm(tmp)
    function create_hypre_vector()
        x = HYPREVector(comm, 1, n)
        asm = HYPRE.start_assemble!(x)
        HYPRE.finish_assemble!(asm)
        return x
    end
    A = HYPREMatrix(comm, 1, n)
    r = create_hypre_vector()
    p = create_hypre_vector()
    p_prec.data[:hypre_system] = (A, r, p)
    return (A, r, p)
end

function JutulDarcy.update_p_rhs!(r_p::HYPRE.HYPREVector, y, ncomp, bz, w_p, p_prec, mode)
    helper = p_prec.data[:assembly_helper]
    inner_hypre_p_rhs!(r_p, y, ncomp, bz, w_p, helper, mode)
end

function inner_hypre_p_rhs!(r_p, y, ncomp, bz, w_p, helper, mode)
    R_p = helper.native_zeroed_buffer
    ix = helper.indices
    if mode == :forward
        @inbounds for i in eachindex(R_p)
            v = 0.0
            for b = 1:ncomp
                v += y[(i-1)*bz + b]*w_p[b, i]
            end
            R_p[i] = v
        end
    else
        @inbounds for i in eachindex(R_p)
            R_p[i] = y[(i-1)*bz + 1]*w_p[1, i]
        end
    end

    Jutul.local_hypre_copy!(r_p, R_p, ix)
    @. R_p = 0.0
end

function get_p_buffer(Δp, p_buf)
    if ismissing(p_buf)
        nvalues = Δp.iupper - Δp.ilower + 1
        p_buf = zeros(nvalues)
    end
    return p_buf::Vector{Float64}
end

function JutulDarcy.correct_residual_and_increment_pressure!(y, x, Δp::HYPRE.HYPREVector, bz, buf, A, p_buf = get_p_buffer(Δp, p_buf))
    get_p_buffer(Δp, p_buf)
    copy!(p_buf, Δp)
    JutulDarcy.correct_residual_and_increment_pressure!(y, x, p_buf, bz, buf, A)
end

function JutulDarcy.increment_pressure!(x, Δp::HYPRE.HYPREVector, bz, p_buf = get_p_buffer(Δp, p_buf))
    copy!(p_buf, Δp)
    JutulDarcy.increment_pressure!(x, p_buf, bz)
end
