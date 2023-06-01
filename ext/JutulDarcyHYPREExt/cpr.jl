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
    @assert n == iupper - ilower + 1 "$(n-1) != $ilower -> $iupper"

    assembler = HYPRE.start_assemble!(A_p)

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
    HYPRE.finish_assemble!(assembler)
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

    if J isa SparseMatrixCSC
        return JutulDarcy.create_pressure_system(nothing, J, n)
    else
        A = HYPREMatrix(comm, 1, n)
        r = create_hypre_vector()
        p = create_hypre_vector()
        p_prec.data[:hypre_system] = (A, r, p)
        return (A, r, p)
    end
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


function Jutul.parray_update_preconditioners!(sim, cpr::CPRPreconditioner{<:BoomerAMGPreconditioner, <:Any}, preconditioners, recorder, tmr)
    offset = sim.storage.process_offset
    n = sim.storage.nc_process
    comm = sim.storage.comm

    function create_hypre_vector()
        x = HYPREVector(comm, offset + 1, offset + n)
        asm = HYPRE.start_assemble!(x)
        HYPRE.finish_assemble!(asm)
        return x
    end
    if isnothing(cpr.A_p)
        # @info "SETTING CPR START" n
        cpr.A_p = HYPREMatrix(comm, offset + 1, offset + n)
        cpr.r_p = create_hypre_vector()
        cpr.p = create_hypre_vector()
        cpr.np = n
    end
    A_p = cpr.A_p
    r_p = cpr.r_p
    x_p = cpr.p

    map(sim.storage.simulators, preconditioners) do sim, prec
        tic!(tmr)
        sys = sim.storage.LinearizedSystem
        model = sim.model
        storage = sim.storage
        prec.A_p = A_p
        prec.p = x_p
        prec.r_p = r_p
        prec.np = n

        prec.pressure_precond.data[:hypre_system] = (A_p, r_p, x_p)
        Jutul.update_preconditioner!(prec, sys, model, storage, recorder, sim.executor)
        toc!(tmr, "Update preconditioner")
        prec
    end

    return (cpr, preconditioners)
end

function Jutul.parray_preconditioner_apply!(global_out, main_prec::CPRPreconditioner{<:BoomerAMGPreconditioner, <:Any}, X, preconditioners, simulator, arg...)
    tmr = simulator.storage.global_timer
    global_cell_vector = simulator.storage.distributed_cell_buffer
    global_buf = simulator.storage.distributed_residual_buffer
    tic!(tmr)
    # toc!(tmr, "Apply preconditioner communication")
    tic!(tmr)
    map(local_values(X), preconditioners, ghost_values(X)) do x, prec, x_g
        @. x_g = 0.0
        JutulDarcy.apply_cpr_first_stage!(prec, x, arg...)
        nothing
    end
    toc!(tmr, "Apply CPR stage 1")
    tic!(tmr)
    # The following is an unsafe version of this:
    # copy!(global_cell_vector, main_prec.p)
    map(own_values(global_cell_vector), preconditioners) do ov, prec
        helper = prec.pressure_precond.data[:assembly_helper]
        indices = helper.indices
        indices::Vector{HYPRE.HYPRE_BigInt}
        nvalues = indices[end] - indices[1] + 1
        HYPRE.@check HYPRE.HYPRE_IJVectorGetValues(main_prec.p, nvalues, indices, ov)
    end
    # End unsafe shenanigans

    # consistent!(global_cell_vector) |> wait
    toc!(tmr, "Apply preconditioner communication")
    tic!(tmr)

    map(own_values(global_buf), own_values(global_cell_vector), preconditioners) do dx, dp, prec
        bz = prec.block_size
        for i in eachindex(dp)
            JutulDarcy.set_dp!(dx, bz, dp, i)
        end
        nothing
    end
    n = length(global_buf)
    if main_prec.full_system_correction
        mul_ix = nothing
    else
        mul_ix = 1
    end
    distributed_mul! = setup_distributed_mul!(tmr, simulator.storage.simulators, mul_ix)
    full_op = LinearOperator(Float64, n, n, false, false, distributed_mul!)
    mul!(X, full_op, global_buf, -1.0, true)

    map(local_values(global_out), local_values(X), preconditioners, local_values(global_cell_vector), ghost_values(X)) do y, x, prec, dp, x_g
        @. x_g = 0.0
        apply!(y, prec.system_precond, x, arg...)
        bz = prec.block_size
        JutulDarcy.increment_pressure!(y, dp, bz)
        nothing
    end
    toc!(tmr, "CPR second stage")
    tic!(tmr)
    consistent!(global_out) |> wait
    toc!(tmr, "Apply preconditioner communication")
    global_out
end
