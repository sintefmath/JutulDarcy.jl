function JutulDarcy.update_pressure_system!(A_p::HYPRE.HYPREMatrix, p_prec, A::Jutul.StaticSparsityMatrixCSR, w_p, bz, ctx, executor)
    D = p_prec.data
    if !haskey(D, :assembly_helper)
        (; ilower, iupper) = A_p
        D[:assembly_helper] = Jutul.generate_hypre_assembly_helper(A, executor, ilower, iupper)
    end
    helper = D[:assembly_helper]
    I_buf, J_buf, V_buffers, = D[:assembly_helper]
    update_pressure_system_hypre!(I_buf, J_buf, V_buffers, A_p, A, w_p, executor, helper.n)
end

function update_pressure_system_hypre!(single_buf, longer_buf, V_buffers, A_p, A, w_p, executor, n)
    nzval = SparseArrays.nonzeros(A)
    cols = Jutul.colvals(A)

    @assert length(single_buf) == 1
    (; iupper, ilower) = A_p
    owned_num = iupper - ilower + 1 # 5
    @assert n == iupper - ilower + 1 "$(n-1) != $ilower -> $iupper"
    offset = ilower - 1

    assembler = HYPRE.start_assemble!(A_p)

    T = eltype(A)
    for row in 1:n
        pos_ix = nzrange(A, row)
        k = length(pos_ix)
        I = single_buf
        I[1] = row + offset
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
            J[ki] = cols[ri] + offset
        end
        HYPRE.assemble!(assembler, I, J, V_buf)
    end
    HYPRE.finish_assemble!(assembler)
end

function JutulDarcy.update_p_rhs!(r_p::HYPRE.HYPREVector, y, bz, w_p)
    tmp = zeros(length(y)÷bz)
    for i in eachindex(tmp)
        v = 0.0
        for b = 1:bz
            v += y[(i-1)*bz + b]*w_p[b, i]
        end
        tmp[i] = v
    end
    (; iupper, ilower) = r_p
    @assert length(tmp) == iupper - ilower + 1
    ix = Int32.(ilower:iupper)

    assembler = HYPRE.start_assemble!(r_p)
    HYPRE.assemble!(assembler, ix, tmp)
    HYPRE.finish_assemble!(assembler)
end

function JutulDarcy.correct_residual_for_dp!(y, x, Δp::HYPRE.HYPREVector, bz, buf, A)
    tmp = zeros(length(y)÷bz)
    copy!(tmp, Δp)
    JutulDarcy.correct_residual_for_dp!(y, x, tmp, bz, buf, A)
end

function JutulDarcy.increment_pressure!(x, Δp::HYPRE.HYPREVector, bz)
    tmp = zeros(length(x)÷bz)
    copy!(tmp, Δp)
    JutulDarcy.increment_pressure!(x, tmp, bz)
end
