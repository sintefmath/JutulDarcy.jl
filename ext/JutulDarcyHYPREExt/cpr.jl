function JutulDarcy.update_pressure_system!(A_p::HYPRE.HYPREMatrix, p_prec, A, w_p, bz, ctx, executor)
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
    @assert length(single_buf) == 1
    (; iupper, ilower) = A_p
    @assert n == iupper - ilower + 1 "$(n-1) != $ilower -> $iupper"
    assembler = HYPRE.start_assemble!(A_p)
    assemble_into_hypre_psystem!(A_p, A, assembler, w_p, single_buf, longer_buf, V_buffers, executor, n)
    HYPRE.finish_assemble!(assembler)
end

function assemble_into_hypre_psystem!(A_p::HYPRE.HYPREMatrix, A::Jutul.StaticSparsityMatrixCSR, assembler, w_p, single_buf, longer_buf, V_buffers, executor, n)
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
            next_value = 0.0
            for component in axes(w_p, 1)
                next_value += w_p[component, row]*A_block[component, 1]
            end
            V_buf[ki] = next_value
            J[ki] = Jutul.executor_index_to_global(executor, col, :column)
            num_added += 1
        end
        HYPRE.assemble!(assembler, I, J, V_buf)
    end
end

function assemble_into_hypre_psystem!(A_p::HYPRE.HYPREMatrix, A::SparseMatrixCSC, assembler, w_p, single_buf, longer_buf, V_buffers, executor, n)
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
            next_value = 0.0
            for component in axes(w_p, 1)
                next_value += w_p[component, row]*A_block[component, 1]
            end
            V_buf[ki] = next_value
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

function JutulDarcy.update_p_rhs!(r_p::HYPRE.HYPREVector, y, bz, w_p, p_prec)
    helper = p_prec.data[:assembly_helper]
    inner_hypre_p_rhs!(r_p, y, bz, w_p, helper)
end

function inner_hypre_p_rhs!(r_p, y, bz, w_p, helper)
    R_p = helper.native_zeroed_buffer
    ix = helper.indices

    @inbounds for i in eachindex(R_p)
        v = 0.0
        for b = 1:bz
            v += y[(i-1)*bz + b]*w_p[b, i]
        end
        R_p[i] = v
    end

    Jutul.local_hypre_copy!(r_p, R_p, ix)
    @. R_p = 0.0
end

function JutulDarcy.correct_residual_for_dp!(y, x, Δp::HYPRE.HYPREVector, bz, buf, A)
    nvalues = Δp.iupper - Δp.ilower + 1
    tmp = zeros(nvalues)
    copy!(tmp, Δp)
    JutulDarcy.correct_residual_for_dp!(y, x, tmp, bz, buf, A)
end

function JutulDarcy.increment_pressure!(x, Δp::HYPRE.HYPREVector, bz)
    nvalues = Δp.iupper - Δp.ilower + 1
    tmp = zeros(nvalues)
    copy!(tmp, Δp)
    JutulDarcy.increment_pressure!(x, tmp, bz)
end


function Jutul.parray_update_preconditioners!(sim::Jutul.PArraySimulator, cpr::CPRPreconditioner{<:BoomerAMGPreconditioner, <:Any}, preconditioners, recorder, tmr)
    offset = sim.storage.process_offset
    n = sim.storage.nc_process
    comm = sim.storage.comm
    if sim.storage[:number_of_processes] > 1
        @assert sim.backend isa Jutul.MPI_PArrayBackend "Cannot use HYPRE with emulated multiple processes."
    end

    function create_hypre_vector()
        x = HYPREVector(comm, offset + 1, offset + n)
        asm = HYPRE.start_assemble!(x)
        HYPRE.finish_assemble!(asm)
        return x
    end
    if isnothing(cpr.A_p)
        cpr.A_p = HYPREMatrix(comm, offset + 1, offset + n)
        cpr.r_p = create_hypre_vector()
        cpr.p = create_hypre_vector()
        cpr.np = n
    end
    A_p = cpr.A_p
    r_p = cpr.r_p
    x_p = cpr.p

    map(sim.storage.simulators, preconditioners) do sim, prec
        # tic!(tmr)
        sys = sim.storage.LinearizedSystem
        model = sim.model
        storage = sim.storage
        prec.A_p = A_p
        prec.p = x_p
        prec.r_p = r_p
        prec.np = n

        prec.pressure_precond.data[:hypre_system] = (A_p, r_p, x_p)
        Jutul.update_preconditioner!(prec, sys, model, storage, recorder, sim.executor)
        prec
    end

    return (cpr, preconditioners)
end
