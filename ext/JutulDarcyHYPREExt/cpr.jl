function JutulDarcy.update_pressure_system!(A_p::HYPRE.HYPREMatrix, p_prec, A::Jutul.StaticSparsityMatrixCSR, w_p, bz, ctx, executor)
    D = p_prec.data
    if !haskey(D, :assembly_helper)
        D[:assembly_helper] = Jutul.generate_hypre_assembly_helper(A, executor)
    end
    helper = D[:assembly_helper]
    I_buf, J_buf, V_buffers, = D[:assembly_helper]
    update_pressure_system_hypre!(I_buf, J_buf, V_buffers, A_p, A, w_p, executor, helper.n)
end

function update_pressure_system_hypre!(single_buf, longer_buf, V_buffers, A_p, A, w_p, executor, n)
    nzval = SparseArrays.nonzeros(A)
    cols = Jutul.colvals(A)

    @assert length(single_buf) == 1
    assembler = HYPRE.start_assemble!(A_p)

    T = eltype(A)
    for row in 1:n
        pos_ix = nzrange(A, row)
        k = length(pos_ix)
        I = single_buf
        I[1] = Jutul.executor_index_to_global(executor, row, :row)
        J = longer_buf
        resize!(J, k)
        V_buf = V_buffers[k]
        @inbounds for ki in 1:k
            ri = pos_ix[ki]
            A_block = nzval[ri]
            next_value = 0.0
            for component in axes(w_p, 1)
                next_value += w_p[component, row]*A_block[component, 1]
            end
            V_buf[ki] = next_value
            J[ki] = Jutul.executor_index_to_global(executor, cols[ri], :column)
        end
        HYPRE.assemble!(assembler, I, J, V_buf)
    end
    HYPRE.finish_assemble!(assembler)
end