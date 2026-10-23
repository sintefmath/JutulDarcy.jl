function hypre_assembly_helper!(p_prec, A, executor, ilower, iupper)
    return get!(p_prec.data, :assembly_helper) do
        Jutul.generate_hypre_assembly_helper(
            A, executor, ilower, iupper)
    end
end

function cprw_hypre_assembly_helper!(p_prec, model, lsys, executor, counts, comm, n)
    A = JutulDarcy.reservoir_jacobian(lsys)
    A isa Jutul.StaticSparsityMatrixCSR || throw(ArgumentError(
        "MPI CPRW requires the CSR backend (backend = :csr)."))
    bz = number_of_components(reservoir_model(model).system)
    well_map = JutulDarcy.cpr_construct_well_reservoir_map(model, lsys, bz)
    pressure_matrix = isnothing(well_map) ? A :
        JutulDarcy.pressure_matrix_from_global_jacobian(A, Float64, lsys, well_map)
    ncell = size(A, 1)
    nwell = size(pressure_matrix, 1) - ncell
    # Wells belong to a single rank, but each rank can own a different number
    # of them. HYPRE needs contiguous ownership including these extra rows.
    pressure_counts = MPI.Allgather(n + nwell, comm)
    pressure_offsets = cumsum([0; pressure_counts])
    cell_offsets = cumsum([0; counts])
    rank = MPI.Comm_rank(comm) + 1
    offset = pressure_offsets[rank]
    npressure = pressure_counts[rank]
    global_indices = Vector{HYPRE.HYPRE_BigInt}(undef, ncell + nwell)
    for cell in 1:ncell
        global_cell = Jutul.executor_index_to_global(executor, cell, :column)
        owner = searchsortedlast(cell_offsets, global_cell - 1)
        global_indices[cell] = pressure_offsets[owner] + global_cell - cell_offsets[owner]
    end
    for well in 1:nwell
        global_indices[ncell + well] = offset + n + well
    end
    # Ghost reservoir rows are skipped; local well rows follow all local cells
    # in the CSR staging matrix and follow owned cells in the HYPRE matrix.
    owned_rows = [1:n; (ncell + 1):(ncell + nwell)]
    positions = reduce(vcat, (collect(nzrange(pressure_matrix, row)) for row in owned_rows); init = Int[])
    columns = Jutul.colvals(pressure_matrix)
    indices = HYPRE.HYPRE_BigInt.((offset + 1):(offset + npressure))
    csr = (
        nrows = HYPRE.HYPRE_Int(npressure),
        ncols = HYPRE.HYPRE_Int[length(nzrange(pressure_matrix, row)) for row in owned_rows],
        rows = indices,
        cols = global_indices[columns[positions]],
        values = zeros(HYPRE.HYPRE_Complex, length(positions)),
        positions = positions,
    )
    p_prec.data[:assembly_helper] = (
        I = Int[], J = Int[], V = nothing,
        indices = indices, native_zeroed_buffer = zeros(npressure),
        n = npressure, csr = csr, nreservoir = n,
        pressure_matrix = pressure_matrix, well_reservoir_map = well_map,
    )
    return offset, npressure
end

function JutulDarcy.update_cpr_pressure_system!(A_p::HYPRE.HYPREMatrix, p_prec,
        A, w_p, ctx, executor, well_map::NamedTuple, model, storage)
    helper = p_prec.data[:assembly_helper]
    pressure_matrix = helper.pressure_matrix
    # Apply the existing CSR reduction and controls before assembling owned rows.
    JutulDarcy.update_cpr_pressure_system!(pressure_matrix, p_prec,
        A, w_p, ctx, executor, well_map, model, storage)
    values = nonzeros(pressure_matrix)
    csr = helper.csr
    for i in eachindex(csr.positions)
        csr.values[i] = values[csr.positions[i]]
    end
    return assemble_pressure_system_hypre_csr!(A_p, csr)
end

function JutulDarcy.update_pressure_system!(A_p::HYPRE.HYPREMatrix, p_prec, A, w_p, ctx, executor, well_reservoir_map::Nothing)
    (; ilower, iupper) = A_p
    helper = hypre_assembly_helper!(
        p_prec, A, executor, ilower, iupper)
    I_buf, J_buf, V_buffers, = helper
    is_adjoint = Val(Jutul.represented_as_adjoint(matrix_layout(ctx)))
    ncomp = Val(size(w_p, 1))
    if A isa Jutul.StaticSparsityMatrixCSR
        update_pressure_system_hypre_csr!(
            A_p, A, w_p, helper, is_adjoint, ncomp)
    else
        update_pressure_system_hypre!(I_buf, J_buf, V_buffers,
            A_p, A, w_p, executor, helper.n, is_adjoint, ncomp)
    end
end

function update_pressure_system_hypre_csr!(A_p, A, w_p, helper,
        ::Val{is_adjoint}, ::Val{ncomp}) where {is_adjoint, ncomp}
    csr = helper.csr
    @assert !isnothing(csr)
    # Only owned rows go to HYPRE; under MPI the local Jacobian also holds ghost rows.
    (; iupper, ilower) = A_p
    n = iupper - ilower + 1
    @assert size(A, 1) >= n
    nzval = SparseArrays.nonzeros(A)
    @inbounds for row in 1:n
        for position in nzrange(A, row)
            csr.values[position] = JutulDarcy.reduce_to_pressure(
                nzval[position], w_p, row, ncomp, is_adjoint)
        end
    end
    return assemble_pressure_system_hypre_csr!(A_p, csr)
end

function assemble_pressure_system_hypre_csr!(A_p, csr)
    assembler = HYPRE.start_assemble!(A_p)
    HYPRE.@check HYPRE.HYPRE_IJMatrixSetValues(
        A_p, csr.nrows, csr.ncols, csr.rows, csr.cols, csr.values)
    HYPRE.finish_assemble!(assembler)
    return A_p
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
    comm = MPI.COMM_SELF
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

function JutulDarcy.update_p_rhs!(r_p::HYPRE.HYPREVector, y, ncomp, bz,
        w_p, pressure_matrix::HYPRE.HYPREMatrix, mode,
        pressure_preconditioner::BoomerAMGPreconditioner)
    helper = pressure_preconditioner.data[:assembly_helper]
    inner_hypre_p_rhs!(r_p, y, ncomp, bz, w_p, helper, mode)
end

function inner_hypre_p_rhs!(r_p, y, ncomp, bz, w_p, helper, mode)
    R_p = helper.native_zeroed_buffer
    ix = helper.indices
    n = get(helper, :nreservoir, length(R_p))
    if mode == :forward
        @inbounds for i in 1:n
            v = 0.0
            for b = 1:ncomp
                v += y[(i-1)*bz + b]*w_p[b, i]
            end
            R_p[i] = v
        end
    else
        @inbounds for i in 1:n
            R_p[i] = y[(i-1)*bz + 1]*w_p[1, i]
        end
    end
    # Auxiliary well entries stay zero: the buffer is cleared after each copy.
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
