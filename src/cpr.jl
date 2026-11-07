
struct CPRStorage{P, R, S, F, W, V, P_v, M}
    "pressure system"
    A_p::P
    "pressure residual"
    r_p::R
    "last pressure approximation"
    p::R
    "buffer of size equal to full system rhs"
    x_ps::S
    "buffer of size equal to full system rhs"
    r_ps::S
    "full system (operator)"
    A_ps::F
    "pressure weights"
    w_p::W
    "weight for right hand side"
    w_rhs::V
    np::Int
    block_size::Int
    number_of_components::Int
    "id for keeping track of usage, typically `objectid` of reservoir jacobian"
    id::UInt64
    "optional pressure-sized buffer"
    p_buffer::P_v
    "optional well+reservoir mappings for CPRW"
    well_reservoir_map::M
    float_type::DataType
end

function allocate_cpr_array(prototype, ::Type{T}, dimensions...) where T
    values = similar(prototype, T, dimensions...)
    fill!(values, zero(T))
    return values
end

function CPRStorage(p_prec, lin_op, full_jac, ncomp = missing;
        T = missing,
        p_buffer::Bool = true,
        well_reservoir_map = nothing,
        lsys = missing
    )
    T_b = eltype(full_jac)
    @assert T_b<:StaticMatrix
    if ismissing(T)
        T = eltype(T_b)
    end
    # TODO: Update the sizes here from well_map if it exists.
    @assert T<:Real
    ncell = size(full_jac, 1)
    if isnothing(well_reservoir_map)
        np = ncell
    else
        np = well_reservoir_map.np
    end
    bz = size(T_b, 1)
    if ismissing(ncomp)
        ncomp = bz
    end
    A_p, r_p, p = create_pressure_system(p_prec, full_jac, lsys, np, T, well_reservoir_map)
    prototype = nonzeros(full_jac)
    solution = allocate_cpr_array(prototype, T, ncell*bz)
    residual = allocate_cpr_array(prototype, T, ncell*bz)
    w_p = allocate_cpr_array(prototype, T, ncomp, np)
    w_rhs = zeros(T, ncomp)
    w_rhs[1] = one(T)
    w_rhs = SVector{ncomp, T}(w_rhs)
    if p_buffer
        p_buf = allocate_cpr_array(prototype, T, np)
    else
        p_buf = missing
    end
    return CPRStorage(A_p, r_p, p, solution, residual, lin_op, w_p, w_rhs, np, bz, ncomp, objectid(full_jac), p_buf, well_reservoir_map, T)
end

function CPRStorage(np::Int, bz::Int, lin_op, psys::Tuple, solution, residual, T = Float64, id = zero(UInt64); ncomp = bz)
    A_p, r_p, p = psys
    w_p = zeros(T, ncomp, np)
    w_rhs = zeros(T, ncomp)
    w_rhs[1] = one(T)
    w_rhs = SVector{ncomp, T}(w_rhs)
    return CPRStorage(A_p, r_p, p, solution, residual, lin_op, w_p, w_rhs, np, bz, ncomp, id, zeros(T, np), nothing, T)
end


mutable struct CPRAMGRebuildStorage{V, T}
    values::V
    squared_norm::T
end

"""
    CPRPreconditioner(p = reservoir_system_amg(), s = ILUZeroPreconditioner(); strategy = :quasi_impes, weight_scaling = :unit, update_frequency = 1, update_interval = :iteration, partial_update = true)

Construct a constrained pressure residual (CPR) preconditioner.

By default, this is a AMG-BILU(0) version (algebraic multigrid for pressure,
block-ILU(0) for the global system).

The default `amg_rebuild_tolerance = 0.05` gates scheduled full updates using
`norm(nonzeros(A_p) - previous_full_values)/norm(previous_full_values)`.
Changes at or below the tolerance take the partial-update path instead,
regardless of the partial-update interval. The reference is saved after each
successful full update; initialization and storage replacement always require
a full update. Set the tolerance to `0.0` to always perform scheduled full
updates without checking matrix changes. This assumes the pressure matrix's
sparsity pattern is fixed.
"""
mutable struct CPRPreconditioner{P, S} <: JutulPreconditioner
    storage::Union{CPRStorage, Nothing}
    pressure_precond::P
    system_precond::S
    strategy::Symbol
    variant::Symbol
    weight_scaling::Symbol
    update_frequency::Int             # Frequency for full pressure-preconditioner updates
    update_interval::Symbol           # iteration, ministep, step, ...
    update_frequency_partial::Int     # Frequency for partial pressure-preconditioner updates
    update_interval_partial::Symbol   # iteration, ministep, step, ...
    partial_update::Bool              # Enable partial pressure-preconditioner updates
    amg_rebuild_tolerance::Float64
    amg_rebuild_storage::Union{Nothing, CPRAMGRebuildStorage}
    full_system_correction::Bool
    p_rtol::Union{Float64, Nothing}
    npre::Int
    npost::Int
    mode::Symbol
    psolver
end

function CPRPreconditioner(p = reservoir_system_amg(), s = ILUZeroPreconditioner();
        strategy = :true_impes,
        variant = :cpr,
        weight_scaling = :unit,
        npre = 1,
        npost = 0,
        update_frequency = 1,
        update_interval = :iteration,
        update_frequency_partial = 1,
        update_interval_partial = :iteration,
        p_rtol = nothing,
        full_system_correction = true,
        partial_update = update_interval == :once,
        amg_rebuild_tolerance = 0.05,
        mode = :forward
    )
    if !(amg_rebuild_tolerance >= 0)
        throw(ArgumentError("amg_rebuild_tolerance must be nonnegative"))
    end
    return CPRPreconditioner(
        nothing,
        p,
        s,
        strategy,
        variant,
        weight_scaling,
        update_frequency,
        update_interval,
        update_frequency_partial,
        update_interval_partial,
        partial_update,
        amg_rebuild_tolerance,
        nothing,
        full_system_correction,
        p_rtol,
        npre,
        npost,
        mode,
        nothing
    )
end

function update_preconditioner!(cpr::CPRPreconditioner, lsys::Jutul.JutulLinearSystem, ctx_outer, model, storage, recorder, executor; update_system_precond = true, T = eltype(Jutul.vector_residual(lsys)))
    rmodel = reservoir_model(model)
    ctx = rmodel.context
    update_p, update_p_partial = update_cpr_internals!(
        cpr, lsys, model, storage, recorder, executor, T)
    # The full-system smoother follows the current Jacobian on every call,
    # independently of how much of the pressure AMG hierarchy is reused.
    if update_system_precond
        @tic "s-precond" update_preconditioner!(cpr.system_precond, lsys, ctx, model, storage, recorder, executor)
    end
    if update_p
        @tic "p-precond" update_preconditioner!(cpr.pressure_precond, cpr.storage.A_p, cpr.storage.r_p, ctx, executor)
        update_cpr_amg_rebuild_reference!(cpr)
    elseif update_p_partial
        @tic "p-precond (partial)" partial_update_preconditioner!(cpr.pressure_precond, cpr.storage.A_p, cpr.storage.r_p, ctx, executor)
    end
end

function initialize_cpr_storage!(cpr, model, lsys, bz, T)
    J = reservoir_jacobian(lsys)
    do_setup = isnothing(cpr.storage) || cpr.storage.id != objectid(J)
    if do_setup
        if cpr.variant == :cprw
            well_reservoir_map = cpr_construct_well_reservoir_map(model, lsys, bz)
        else
            well_reservoir_map = nothing
        end
        if cpr.full_system_correction
            op = linear_operator(lsys)
        else
            op = linear_operator(lsys[1,1])
        end
        cpr.storage = CPRStorage(cpr.pressure_precond, op, J, bz,
            well_reservoir_map = well_reservoir_map,
            lsys = lsys,
            T = T
        )
        cpr.amg_rebuild_storage = nothing
    end
    return do_setup
end

# GPUArrays fuses multi-array mapreduce into its existing reduction kernels.
function cpr_amg_rebuild_difference(values, reference)
    squared_difference(a, b) = abs2(a - b)
    initial_sum = zero(eltype(values))
    difference_squared = mapreduce(squared_difference, +, values, reference; init = initial_sum)
    return difference_squared
end

# Base's multi-array mapreduce materializes map(f, a, b), so sum over indices
# for ordinary CPU vectors instead.
function cpr_amg_rebuild_difference(values::Vector, reference::Vector)
    indices = eachindex(values, reference)
    initial_sum = zero(eltype(values))
    difference_squared = sum(indices; init = initial_sum) do i
        @inbounds difference = values[i] - reference[i]
        return abs2(difference)
    end
    return difference_squared
end

function cpr_amg_rebuild_change(values, reference::CPRAMGRebuildStorage)
    reference_values = reference.values
    difference_squared = cpr_amg_rebuild_difference(values, reference_values)
    reference_squared = reference.squared_norm
    finite_norms = isfinite(difference_squared) && isfinite(reference_squared)
    if !finite_norms
        return Inf
    elseif iszero(reference_squared)
        return iszero(difference_squared) ? 0.0 : Inf
    else
        relative_change_squared = difference_squared/reference_squared
        return sqrt(relative_change_squared)
    end
end

function cpr_amg_rebuild_compatible(values, reference::CPRAMGRebuildStorage)
    reference_values = reference.values
    same_type = typeof(values) === typeof(reference_values)
    same_axes = axes(values) == axes(reference_values)
    return same_type && same_axes
end

function cpr_amg_rebuild_supported(cpr)
    is_csr = cpr.storage.A_p isa Jutul.StaticSparsityMatrixCSR
    is_hypre = cpr.pressure_precond isa BoomerAMGPreconditioner
    return is_csr && !is_hypre
end

function cpr_amg_update_decision(cpr, full, partial)
    tolerance = cpr.amg_rebuild_tolerance
    check_change = full && tolerance > 0 && cpr_amg_rebuild_supported(cpr)
    if check_change
        reference = cpr.amg_rebuild_storage
        pressure_matrix = cpr.storage.A_p
        values = nonzeros(pressure_matrix)
        compatible_reference = !isnothing(reference) && cpr_amg_rebuild_compatible(values, reference)
        if compatible_reference
            change = cpr_amg_rebuild_change(values, reference)
            within_tolerance = isfinite(change) && change <= tolerance
            if within_tolerance
                full = false
                partial = true
            end
        end
    end
    return full, partial
end

function update_cpr_amg_rebuild_reference!(cpr)
    tolerance = cpr.amg_rebuild_tolerance
    if iszero(tolerance) || !cpr_amg_rebuild_supported(cpr)
        cpr.amg_rebuild_storage = nothing
        return nothing
    end
    pressure_matrix = cpr.storage.A_p
    values = nonzeros(pressure_matrix)
    reference = cpr.amg_rebuild_storage
    needs_new_reference = isnothing(reference) || !cpr_amg_rebuild_compatible(values, reference)
    if needs_new_reference
        reference_values = copy(values)
        reference_squared = sum(abs2, values)
        cpr.amg_rebuild_storage = CPRAMGRebuildStorage(reference_values, reference_squared)
    else
        reference_values = reference.values
        copyto!(reference_values, values)
        reference_squared = sum(abs2, values)
        reference.squared_norm = reference_squared
    end
    return nothing
end

function pressure_matrix_from_global_jacobian(J::SparseMatrixCSC, T, lsys, well_reservoir_map::Nothing)
    nzval = zeros(T, nnz(J))
    n = size(J, 2)
    return SparseMatrixCSC(n, n, J.colptr, J.rowval, nzval)
end

function pressure_matrix_from_global_jacobian(J::Jutul.StaticSparsityMatrixCSR, T, lsys, well_reservoir_map::Nothing)
    nzval = allocate_cpr_array(nonzeros(J), T, nnz(J))
    n = size(J, 2)
    # Assume symmetry in sparse pattern, but not values.
    return Jutul.StaticSparsityMatrixCSR(
        nzval, Jutul.colvals(J), J.rowptr, n, n, J.backend;
        nthreads = J.nthreads, minbatch = J.minbatch,
        thread_type = J.thread_type
    )
end

function pressure_matrix_from_global_jacobian(sys_jac::Jutul.StaticSparsityMatrixCSR, T, lsys, well_reservoir_map)
    n = well_reservoir_map.np
    I = well_reservoir_map.I
    J = well_reservoir_map.J
    @assert length(I) == length(J)
    V = zeros(T, length(I))
    A = sparse(J, I, V, n, n)
    A_p_host = Jutul.StaticSparsityMatrixCSR(A;
        nthreads = sys_jac.nthreads,
        minbatch = sys_jac.minbatch
    )
    # Align reservoir variables
    nzmap_reservoir = well_reservoir_map.nzmap_11
    resize!(nzmap_reservoir, nnz(sys_jac))
    ncell = size(sys_jac, 1)
    nwell = size(A_p_host, 1) - ncell
    system_rowptr = Array(sys_jac.rowptr)
    system_colval = Array(Jutul.colvals(sys_jac))
    for row in 1:ncell
        for j in system_rowptr[row]:(system_rowptr[row + 1] - 1)
            col = system_colval[j]
            nzmap_reservoir[j] = find_sparse_position(A_p_host, row, col)
        end
    end
    # Well to well
    nzmap_well = well_reservoir_map.nzmap_22
    for w in 1:nwell
        ix = find_sparse_position(A_p_host, ncell + w, ncell + w)
        push!(nzmap_well, ix)
    end

    # Align well to reservoir variables
    nzmap_12 = well_reservoir_map.nzmap_12
    nzmap_21 = well_reservoir_map.nzmap_21
    for w in 1:nwell
        for cell in well_reservoir_map.well_cells[w]
            ix_12 = find_sparse_position(A_p_host, cell, ncell + w)
            push!(nzmap_12, ix_12)
            ix_21 = find_sparse_position(A_p_host, ncell + w, cell)
            push!(nzmap_21, ix_21)
        end
    end
    A_p = Jutul.transfer_csr_to_backend(sys_jac, A_p_host)
    prepare_cprw_backend_maps!(well_reservoir_map, sys_jac)
    return A_p
end

function cpr_backend_copy(prototype, source::AbstractVector{T}) where T
    if prototype isa Vector
        return source
    end
    destination = similar(prototype, T, length(source))
    copyto!(destination, source)
    return destination
end

function prepare_cprw_backend_maps!(map, system_matrix)
    flatten(groups) = reduce(vcat, groups; init = eltype(groups)())
    map_12 = flatten(map.map_12)
    map_21 = flatten(map.map_21)
    map_22 = flatten(map.map_22)
    if isempty(map.well_cells)
        cells_12 = Int[]
    else
        cells_12 = reduce(vcat,
            (map.well_cells[index] for index in eachindex(map.map_12));
            init = Int[])
    end
    reservoir_cells = size(system_matrix, 1)
    wells_21 = reduce(vcat,
        (fill(index + reservoir_cells, length(map.map_21[index]))
            for index in eachindex(map.map_21));
        init = Int[])
    wells_22 = reduce(vcat,
        (fill(index + reservoir_cells, length(map.map_22[index]))
            for index in eachindex(map.map_22)); init = Int[])
    positions_22 = reduce(vcat,
        (fill(map.nzmap_22[index], length(map.map_22[index]))
            for index in eachindex(map.map_22)); init = Int[])
    prototype = system_matrix.rowptr
    map.backend_data[] = (
        cross_12 = cprw_backend_reduction(prototype, map.nzmap_12, map_12, cells_12),
        cross_21 = cprw_backend_reduction(prototype, map.nzmap_21, map_21, wells_21),
        cross_22 = cprw_backend_reduction(prototype, positions_22, map_22, wells_22),
        cells_12 = cpr_backend_copy(prototype, cells_12),
        well_cell_ptr = cpr_backend_copy(prototype, cumsum([1; length.(map.well_cells)])),
        nzmap_11 = cpr_backend_copy(prototype, map.nzmap_11),
    )
    return map
end

function cprw_backend_reduction(prototype, positions, blocks, entities)
    # One thread sums all fine blocks contributing to a coarse entry. This
    # includes segment diagonals and multiple segments perforating one cell.
    order = sortperm(positions)
    sorted = positions[order]
    starts = findall(i -> i == 1 || sorted[i] != sorted[i - 1], eachindex(sorted))
    return (
        positions = cpr_backend_copy(prototype, sorted[starts]),
        offsets = cpr_backend_copy(prototype, [starts; length(sorted) + 1]),
        blocks = cpr_backend_copy(prototype, blocks[order]),
        entities = cpr_backend_copy(prototype, entities[order]),
    )
end

function create_pressure_system(p_prec, J, lsys, n, T, well_reservoir_map)
    J_p = pressure_matrix_from_global_jacobian(J, T, lsys, well_reservoir_map)
    prototype = nonzeros(J_p)
    return (
        J_p,
        allocate_cpr_array(prototype, T, n),
        allocate_cpr_array(prototype, T, n)
    )
end

function update_cpr_internals!(cpr::CPRPreconditioner, lsys, model, storage, recorder, executor, T)
    A = reservoir_jacobian(lsys)
    rmodel = reservoir_model(model)
    bz = number_of_components(rmodel.system)
    storage_changed = initialize_cpr_storage!(cpr, model, lsys, bz, T)
    do_p_update = storage_changed || should_update_cpr(cpr, recorder, :amg)
    do_p_update_partial = !do_p_update && should_update_cpr(cpr, recorder, :partial)
    ps = rmodel.primary_variables[:Pressure].scale
    # The pressure operator and weights always follow the current Jacobian.
    # The update intervals only select how its preconditioner is refreshed.
    rmodel = reservoir_model(model)
    ctx = rmodel.context
    @tic "weights" w_p = update_weights!(cpr, cpr.storage, model, storage, A, ps)
    cpr_s = cpr.storage
    A_p = cpr_s.A_p
    w_p = cpr_s.w_p
    rw_map = cpr_s.well_reservoir_map
    @tic "pressure system" update_cpr_pressure_system!(A_p, cpr.pressure_precond, A, w_p, ctx, executor, rw_map, model, storage)
    return cpr_amg_update_decision(cpr, do_p_update, do_p_update_partial)
end

function update_cpr_pressure_system!(A_p, p_prec, A, w_p, ctx, executor, rw_map, model, storage)
    if !isnothing(rw_map)
        refresh_cprw_source_values!(rw_map)
    end
    update_pressure_system!(A_p, p_prec, A, w_p, ctx, executor, rw_map)
    if !isnothing(rw_map)
        cprw_fix_bhp_controlled_wells!(A_p, rw_map, model, storage)
    end
    return A_p
end

# A BHP-controlled well's pressure is fixed by its control equation, which is not part of the
# CPRW pressure system. Its mass-balance row would make it a free, nearly volume-less node, so
# keep only the diagonal (the row's rhs is already zero): the pressure stage then keeps BHP fixed.
function cprw_fix_bhp_controlled_wells!(A_p::Jutul.StaticSparsityMatrixCSR, rw_map, model, storage)
    model isa MultiModel || return A_p
    offset = rw_map.np - sum(length, rw_map.well_ranges; init = 0)  # wells come last
    for (key, range) in zip(rw_map.wells, rw_map.well_ranges)
        fkey = haskey(model.models, :Facility) ? :Facility : Symbol(key, :_ctrl)
        haskey(storage, fkey) || continue
        # A host-evaluated facility owns the current controls; its backend mirror
        # only receives numerical state updates.
        facility_storage, = Jutul.submodel_evaluation_pair(storage, model, fkey)
        cfg = facility_storage.state.WellGroupConfiguration
        well = physical_representation(model.models[key].domain)
        first_row = offset + first(range)
        info = well.multiwell
        multisegment = well isa MultiSegmentWell
        if isnothing(info)
            cprw_fix_bhp_controlled_well!(A_p, cfg, key, first_row, length(range), multisegment)
        else
            # Merged wells retain their original control names and node ranges.
            for (index, name) in enumerate(info.names)
                start = multisegment ? index : info.nodes.pos[index]
                count = multisegment ? 1 : info.nodes.pos[index + 1] - start
                cprw_fix_bhp_controlled_well!(A_p, cfg, name, first_row + start - 1, count, multisegment)
            end
        end
    end
    return A_p
end

function cprw_fix_bhp_controlled_well!(A_p, cfg, name, first_row, count, unit_diagonal = false)
    haskey(cfg.operating_controls, name) || return A_p
    if operating_control(cfg, name).target isa BottomHolePressureTarget
        cprw_zero_offdiagonals!(A_p, first_row, count, unit_diagonal)
    end
    return A_p
end

function cprw_zero_offdiagonals!(A_p, first_row, count, unit_diagonal = false)
    nz = nonzeros(A_p)
    cols = Jutul.colvals(A_p)
    for row in first_row:(first_row + count - 1)
        for pos in nzrange(A_p, row)
            if cols[pos] != row
                nz[pos] = zero(eltype(nz))
            elseif unit_diagonal
                nz[pos] = one(eltype(nz))
            end
        end
    end
    return A_p
end

cprw_fix_bhp_controlled_wells!(A_p, rw_map, model, storage) = A_p

# CSC (default) version
function update_pressure_system!(A_p, p_prec, A, w_p, ctx, executor, ::Nothing)
    nz = nonzeros(A_p)
    nz_s = nonzeros(A)
    rows = rowvals(A_p)
    @assert size(nz) == size(nz_s)
    n = A.n
    # Update the pressure system with the same pattern in-place
    tb = minbatch(ctx, n)
    ncomp = size(w_p, 1)
    N = Val(ncomp)
    is_adjoint = Val(Jutul.represented_as_adjoint(matrix_layout(ctx)))
    @batch minbatch=tb for col in 1:n
        update_row_csc!(nz, A_p, w_p, rows, nz_s, col, N, is_adjoint)
    end
end

function update_row_csc!(nz, A_p, w_p, rows, nz_s, col, ::Val{Ncomp}, ::Val{adjoint}) where {Ncomp, adjoint}
    @inbounds for j in nzrange(A_p, col)
        row = rows[j]
        Ji = nz_s[j]
        nz[j] = reduce_to_pressure(Ji, w_p, col, Ncomp, adjoint)
    end
end

# CSR version
function update_pressure_system!(A_p::Jutul.StaticSparsityMatrixCSR, p_prec, A::Jutul.StaticSparsityMatrixCSR, w_p, ctx, executor, ::Nothing)
    T_p = eltype(A_p)
    nz = nonzeros(A_p)
    nz_s = nonzeros(A)
    cols = Jutul.colvals(A)
    @assert size(nz) == size(nz_s)
    n = size(A_p, 1)
    # Update the pressure system with the same pattern in-place
    tb = minbatch(ctx, n)
    ncomp = size(w_p, 1)
    N = Val(ncomp)
    is_adjoint = Val(Jutul.represented_as_adjoint(matrix_layout(ctx)))
    update_rows_csr(nz, A_p, w_p, cols, nz_s, N, is_adjoint, n, tb)
end

function update_rows_csr(nz, A_p, w_p, cols, nz_s, N, is_adjoint, n, tb)
    @batch minbatch=tb for row in 1:n
        update_row_csr!(nz, A_p, w_p, cols, nz_s, row, N, is_adjoint)
    end
end

function update_row_csr!(nz, A_p, w_p, cols, nz_s, row, ::Val{Ncomp}, ::Val{adjoint}) where {Ncomp, adjoint}
    @inbounds for j in nzrange(A_p, row)
        Ji = nz_s[j]
        nz[j] = reduce_to_pressure(Ji, w_p, row, Ncomp, adjoint)
    end
end

function update_pressure_system!(A_p::Jutul.StaticSparsityMatrixCSR, p_prec, A_s::Jutul.StaticSparsityMatrixCSR, w_p, ctx, executor, well_reservoir_map)
    # CPRW has a different sparsity pattern than the reservoir matrix, so we
    # have to do some extra work to get everything in the right spots.
    T_p = eltype(A_p)
    nz_p = nonzeros(A_p)
    nz_s = nonzeros(A_s)
    cols = Jutul.colvals(A_s)
    np = size(A_p, 1)
    ncells = size(A_s, 1)
    @assert np == well_reservoir_map.np
    tb = minbatch(ctx, np)
    ncomp = size(w_p, 1)
    N = Val(ncomp)
    is_adjoint = Val(Jutul.represented_as_adjoint(matrix_layout(ctx)))
    cprw_update_pressure_system!(nz_p, nz_s, A_s, A_p, w_p, well_reservoir_map, tb, ncells, np, N, is_adjoint)
    return A_p
end

function cprw_update_pressure_system!(nz_p, nz_s, A_s, A_p, w_p, well_reservoir_map, tb, ncells, np,::Val{ncomp}, ::Val{is_adjoint}) where {ncomp, is_adjoint}
    (; map_12, map_21, map_22, nzmap_11, nzmap_12, nzmap_21, nzmap_22) = well_reservoir_map
    for positions in (nzmap_12, nzmap_21, nzmap_22), position in positions
        nz_p[position] = zero(eltype(nz_p))
    end
    # Reservoir internal connections
    @assert length(nzmap_11) == length(nz_s)
    @batch minbatch=tb for cell in 1:ncells
        @inbounds for pos_s in nzrange(A_s, cell)
            pos_p = nzmap_11[pos_s]
            Ji = nz_s[pos_s]
            nz_p[pos_p] = reduce_to_pressure(Ji, w_p, cell, ncomp, is_adjoint)
        end
    end
    # Reservoir differentiated wrt wells
    nz_12 = well_reservoir_map.nzval_12
    ix = 1
    @inbounds for i in eachindex(map_12)
        map_12_i = map_12[i]
        wc = well_reservoir_map.well_cells[i]
        # nzmap_12_i = nzmap_12[i]
        @inbounds for j in eachindex(map_12_i)
            cell = wc[j]
            Ji = cprw_block(nz_12, map_12_i[j])
            nz_p[nzmap_12[ix]] += reduce_to_pressure(Ji, w_p, cell, ncomp, is_adjoint)
            ix += 1
        end
    end
    # Wells differentiated wrt reservoir
    nz_21 = well_reservoir_map.nzval_21
    ix = 1
    @inbounds for i in eachindex(map_21)
        map_21_i = map_21[i]
        @inbounds for j in eachindex(map_21_i)
            well = i + ncells
            Ji = cprw_block(nz_21, map_21_i[j])
            nz_p[nzmap_21[ix]] += reduce_to_pressure(Ji, w_p, well, ncomp, is_adjoint)
            ix += 1
        end
    end
    # Wells differentiated wrt wells
    nz_22 = well_reservoir_map.nzval_22
    @inbounds for i in eachindex(map_22, nzmap_22)
        map_22_i = map_22[i]
        @inbounds for j in eachindex(map_22_i)
            well = i + ncells
            Ji = cprw_block(nz_22, map_22_i[j])
            nz_p[nzmap_22[i]] += reduce_to_pressure(Ji, w_p, well, ncomp, is_adjoint)
        end
    end
    return nz_p
end

@inline cprw_block(values, positions) = map(pos -> iszero(pos) ? zero(eltype(values)) : values[pos], positions)

function reduce_to_pressure(Ji, w_p, cell, Ncomp, adjoint)
    out = 0.0
    @inbounds for b = 1:Ncomp
        if adjoint
            out += Ji[1, b]*w_p[b, cell]
        else
            out += Ji[b, 1]*w_p[b, cell]
        end
    end
    return out
end

function operator_nrows(cpr::CPRPreconditioner)
    return size(cpr.storage.A_ps, 1)
end

using Krylov
function apply!(x, cpr::CPRPreconditioner, r0, arg...)
    cpr_s = cpr.storage
    # Get buffers and set working values
    r = cpr_s.r_ps
    @. r  = r0
    buf = cpr_s.x_ps
    A_ps = cpr_s.A_ps
    smoother = cpr.system_precond
    bz = cpr_s.block_size
    # Zero out buffer, just in case (assumed by some solvers)
    fill!(x, zero(eltype(x)))
    # presmooth
    if cpr.npre > 0
        @tic "cpr smoother" apply_cpr_smoother!(x, r, buf, smoother, A_ps, cpr.npre)
    end
    @tic "cpr pressure stage" apply_cpr_pressure_stage!(cpr, cpr_s, r, arg...)
    # postsmooth
    if cpr.npost > 0
        correct_residual_and_increment_pressure!(r, x, cpr_s.p, bz, buf, cpr_s.A_ps, cpr_s.p_buffer)
        @tic "cpr smoother" apply_cpr_smoother!(x, r, buf, smoother, A_ps, cpr.npost, skip_last = true)
    else
        @tic "p increment" increment_pressure!(x, cpr_s.p, bz, cpr_s.p_buffer)
    end
end

function apply_cpr_smoother!(x, r, buf, smoother, A_ps, n; skip_last = false)
    for i in 1:n
        apply!(buf, smoother, r)
        @. x += buf
        if i < n || !skip_last
            correct_residual!(r, A_ps, buf)
        end
    end
end

function correct_residual!(r, A, x)
    @tic "residual correction" mul!(r, A, x,
        -one(eltype(r)), one(eltype(r)))
end

function apply_cpr_pressure_stage!(cpr::CPRPreconditioner, cpr_s::CPRStorage, r, arg...)
    r_p, w_p, bz, Δp = cpr_s.r_p, cpr_s.w_p, cpr_s.block_size, cpr_s.p
    ncomp = cpr_s.number_of_components
    @tic "p rhs" update_p_rhs!(
        r_p, r, ncomp, bz, w_p, cpr_s.A_p, cpr.mode,
        cpr.pressure_precond)
    # Apply preconditioner to pressure part
    @tic "p apply" begin
        p_rtol = cpr.p_rtol
        p_precond = cpr.pressure_precond
        cpr_p_apply!(Δp, cpr, p_precond, r_p, p_rtol)
    end
end

function update_p_rhs!(r_p, y, ncomp, bz, w_p, pressure_matrix, mode,
        pressure_preconditioner)
    return update_p_rhs!(
        r_p, y, ncomp, bz, w_p, pressure_matrix, mode)
end

function cpr_p_apply!(Δp, cpr, p_precond, r_p, p_rtol)
    apply!(Δp, p_precond, r_p)
    if !isnothing(p_rtol)
        A_p = cpr.storage.A_p
        if isnothing(cpr.psolver)
            cpr.psolver = FgmresWorkspace(A_p, r_p)
        end
        psolve = cpr.psolver
        T = eltype(r_p)
        M = Jutul.PrecondWrapper(linear_operator(p_precond, T,
            nothing, nothing, nothing, nothing, nothing))
        fgmres!(psolve, A_p, r_p, Δp, M = M,
            rtol = convert(T, p_rtol), atol = convert(T, 1e-12),
            itmax = 20)
        @. Δp = psolve.x
    end
end

reservoir_residual(lsys) = lsys.r
reservoir_jacobian(lsys) = lsys.jac

function reservoir_residual(lsys::MultiLinearizedSystem)
    return lsys[1, 1].r
end

function reservoir_jacobian(lsys::MultiLinearizedSystem)
    return lsys[1, 1].jac
end

function update_weights!(cpr, cpr_storage::CPRStorage, model, storage, J, ps)
    np = cpr_storage.np
    w = cpr_storage.w_p
    rmodel = reservoir_model(model)
    rstorage = reservoir_storage(model, storage)
    nc = cpr_weights_for_reservoir!(rmodel, J, cpr, cpr_storage, rstorage, 0)
    if model isa MultiModel && haskey(model.models, :Fractures)
        fmodel = model.models[:Fractures]
        fstorage = storage[:Fractures]
        nf = cpr_weights_for_reservoir!(fmodel, J, cpr, cpr_storage, fstorage, nc)
    else
        nf = 0
    end
    last_porous_media_index = nf + nc
    if np > last_porous_media_index
        wr_map = cpr_storage.well_reservoir_map
        ncomp = size(w, 1)
        scaling = cpr.weight_scaling
        @assert !isnothing(wr_map)
        r = cpr_storage.w_rhs
        for (well, well_range) in zip(wr_map.wells, wr_map.well_ranges)
            well_model = model.models[well]
            if physical_representation(well_model.domain) isa MultiSegmentWell
                # Use the same averaged restriction for every segment of a
                # physical well, so internal segment fluxes cancel in its row.
                cprw_average_weights!(w, wr_map, well_range, last_porous_media_index, rmodel.context)
                continue
            end
            first_column = last_porous_media_index + first(well_range)
            last_column = last_porous_media_index + last(well_range)
            w_i = view(w, :, first_column:last_column)
            nwell_cells = length(well_range)
            wstate = storage[well].state
            if cpr.strategy == :analytical
                update_analytical_cpr_weights!(
                    w_i, well_model, wstate, nothing, nwell_cells, ncomp, scaling,
                    rmodel.context)
            else
                acc_i = wstate.TotalMasses
                update_true_impes_weights!(
                    w_i, acc_i, r, nwell_cells, ncomp, ps, scaling,
                    rmodel.context)
            end
        end
    end
    return w
end

function cprw_average_weights!(weights, map, well_range, ncell, context)
    for well in well_range
        cells = map.well_cells[well]
        for component in axes(weights, 1)
            value = zero(eltype(weights))
            for cell in cells
                value += weights[component, cell]
            end
            weights[component, ncell + well] = isempty(cells) ? one(value) : value/length(cells)
        end
    end
    return weights
end

function cpr_weights_for_reservoir!(model::SimulationModel, J, cpr, cpr_storage, storage, offset)
    np = cpr_storage.np
    nc = number_of_cells(model.domain)
    n = min(np, nc)
    r = cpr_storage.w_rhs
    w = view(cpr_storage.w_p, :, (offset+1):(offset+n))
    ncomp = size(w, 1)
    scaling = cpr.weight_scaling
    strategy = cpr.strategy
    if strategy == :true_impes
        eq_s = storage.equations[:mass_conservation]
        if eq_s isa ConservationLawTPFAStorage
            acc = eq_s.accumulation.entries
            ps = model.primary_variables[:Pressure].scale
        else
            acc = storage.state.TotalMasses
            # This term isn't scaled by dt, so use simple weights instead
            ps = 1.0
        end
        update_true_impes_weights!(
            w, acc, r, n, ncomp, ps, scaling, model.context)
    elseif strategy == :analytical
        rstate = storage.state
        update_analytical_cpr_weights!(
            w, model, rstate, r, n, ncomp, scaling, model.context)
    elseif strategy == :quasi_impes
        update_quasi_impes_weights!(
            w, J, r, n, ncomp, scaling, offset, model.context)
    elseif strategy == :none
        # Do nothing. Already set to one.
    else
        error("Unsupported strategy $(strategy)")
    end
    return n
end

function update_true_impes_weights!(
        w, acc, r, n, ncomp, pressure_scale, scaling, context)
    return true_impes!(w, acc, r, n, ncomp, pressure_scale, scaling)
end

function update_quasi_impes_weights!(
        w, J, r, n, ncomp, scaling, offset, context)
    return quasi_impes!(w, J, r, n, ncomp, scaling, offset)
end

function update_analytical_cpr_weights!(
        w, model, state, r, n, ncomp, scaling, context)
    return cpr_weights_no_partials!(w, model, state, r, n, ncomp, scaling)
end

function true_impes!(w, acc, r, n, bz, arg...)
    if bz == 2
        # Hard coded variants
        true_impes_2!(w, acc, r, n, bz, arg...)
    elseif bz == 3
        true_impes_3!(w, acc, r, n, bz, arg...)
    elseif bz == 4
        true_impes_4!(w, acc, r, n, bz, arg...)
    elseif bz == 5
        true_impes_5!(w, acc, r, n, bz, arg...)
    elseif bz == 8
        true_impes_8!(w, acc, r, n, bz, arg...)
    else
        true_impes_gen!(w, acc, r, n, Val(bz), arg...)
    end
end

function true_impes_2!(w, acc, r, n, bz, p_scale, scaling)
    r_p = SVector{2}(r)
    @inbounds for cell in 1:n
        W = acc[1, cell]
        O = acc[2, cell]

        ∂W∂p = W.partials[1]*p_scale
        ∂O∂p = O.partials[1]*p_scale

        ∂W∂s = W.partials[2]
        ∂O∂s = O.partials[2]

        A = @SMatrix [∂W∂p ∂O∂p; 
                      ∂W∂s ∂O∂s]
        invert_w!(w, A, r_p, cell, bz, scaling)
    end
end

function M_entry(acc, i, j, c)
    return @inbounds acc[j, c].partials[i]
end

# TODO: Turn these into @generated functions

function true_impes_3!(w, acc, r, n, bz, s, scaling)
    r_p = SVector{3}(r)
    f(i, j, c) = M_entry(acc, i, j, c)
    for c in 1:n
        A = @SMatrix    [s*f(1, 1, c) s*f(1, 2, c) s*f(1, 3, c);
                           f(2, 1, c) f(2, 2, c) f(2, 3, c);
                           f(3, 1, c) f(3, 2, c) f(3, 3, c)]
        invert_w!(w, A, r_p, c, bz, scaling)
    end
end

function true_impes_4!(w, acc, r, n, bz, s, scaling)
    r_p = SVector{4}(r)
    f(i, j, c) = M_entry(acc, i, j, c)
    for c in 1:n
        A = @SMatrix    [s*f(1, 1, c) s*f(1, 2, c) s*f(1, 3, c) s*f(1, 4, c);
                           f(2, 1, c) f(2, 2, c) f(2, 3, c) f(2, 4, c);
                           f(3, 1, c) f(3, 2, c) f(3, 3, c) f(3, 4, c);
                           f(4, 1, c) f(4, 2, c) f(4, 3, c) f(4, 4, c)]
        invert_w!(w, A, r_p, c, bz, scaling)
    end
end

function true_impes_5!(w, acc, r, n, bz, s, scaling)
    r_p = SVector{5}(r)
    f(i, j, c) = M_entry(acc, i, j, c)
    for c in 1:n
        A = @SMatrix    [s*f(1, 1, c) s*f(1, 2, c) s*f(1, 3, c) s*f(1, 4, c) s*f(1, 5, c);
                         f(2, 1, c) f(2, 2, c) f(2, 3, c) f(2, 4, c) f(2, 5, c);
                         f(3, 1, c) f(3, 2, c) f(3, 3, c) f(3, 4, c) f(3, 5, c);
                         f(4, 1, c) f(4, 2, c) f(4, 3, c) f(4, 4, c) f(4, 5, c);
                         f(5, 1, c) f(5, 2, c) f(5, 3, c) f(5, 4, c) f(5, 5, c)]
        invert_w!(w, A, r_p, c, bz, scaling)
    end
end

function true_impes_8!(w, acc, r, n, bz, s, scaling)
    r_p = SVector{8}(r)
    f(i, j, c) = M_entry(acc, i, j, c)
    for c in 1:n
        A = @SMatrix    [s*f(1, 1, c) s*f(1, 2, c) s*f(1, 3, c) s*f(1, 4, c) s*f(1, 5, c) s*f(1, 6, c) s*f(1, 7, c) s*f(1, 8, c);
                         f(2, 1, c) f(2, 2, c) f(2, 3, c) f(2, 4, c) f(2, 5, c) f(2, 6, c) f(2, 7, c) f(2, 8, c);
                         f(3, 1, c) f(3, 2, c) f(3, 3, c) f(3, 4, c) f(3, 5, c) f(3, 6, c) f(3, 7, c) f(3, 8, c);
                         f(4, 1, c) f(4, 2, c) f(4, 3, c) f(4, 4, c) f(4, 5, c) f(4, 6, c) f(4, 7, c) f(4, 8, c);
                         f(5, 1, c) f(5, 2, c) f(5, 3, c) f(5, 4, c) f(5, 5, c) f(5, 6, c) f(5, 7, c) f(5, 8, c);
                         f(6, 1, c) f(6, 2, c) f(6, 3, c) f(6, 4, c) f(6, 5, c) f(6, 6, c) f(6, 7, c) f(6, 8, c);
                         f(7, 1, c) f(7, 2, c) f(7, 3, c) f(7, 4, c) f(7, 5, c) f(7, 6, c) f(7, 7, c) f(7, 8, c);
                         f(8, 1, c) f(8, 2, c) f(8, 3, c) f(8, 4, c) f(8, 5, c) f(8, 6, c) f(8, 7, c) f(8, 8, c)]
        invert_w!(w, A, r_p, c, bz, scaling)
    end
end

function true_impes_gen!(w, acc, r, n, ::Val{bz}, p_scale, scaling) where bz
    r_p = SVector{bz}(r)
    r_T = eltype(r)
    A = MMatrix{bz, bz, r_T}(undef)
    for cell in 1:n
        @inbounds for i = 1:bz
            v = acc[i, cell]
            @inbounds A[1, i] = v.partials[1]*p_scale
            @inbounds for j = 2:bz
                A[j, i] = v.partials[j]
            end
        end
        invert_w!(w, SMatrix{bz, bz, r_T}(A), r_p, cell, bz, scaling)
    end
end

function quasi_impes!(w, J, r, n, bz, scaling, offset = 0)
    r_p = SVector{bz}(r)
    @batch for cell = 1:n
        J_b = J[cell + offset, cell + offset]'
        invert_w!(w, J_b, r_p, cell, bz, scaling)
    end
end

@inline function invert_w!(w, J, r, cell, bz, scaling)
    tmp = J\r
    if scaling == :unit
        s = 1.0/norm(tmp)
    else
        s = 1.0
    end
    @inbounds for i = 1:bz
        w[i, cell] = tmp[i]*s
    end
end

function update_p_rhs!(r_p, y, ncomp, bz, w_p, p_prec, mode)
    n = length(y)÷bz
    if mode == :forward
        @batch minbatch = 1000 for i in 1:n
            v = 0.0
            @inbounds for b = 1:ncomp
                v += y[(i-1)*bz + b]*w_p[b, i]
            end
            @inbounds r_p[i] = v
        end
    else
        @batch minbatch = 1000 for i in 1:n
            @inbounds r_p[i] = y[(i-1)*bz + 1]*w_p[1, i]
        end
    end
    for i in (n+1):length(r_p)
        r_p[i] = 0.0
    end
end

function correct_residual_and_increment_pressure!(r, x, Δp, bz, buf, A_ps, p_buf = missing)
    # x = x' + Δx
    # A (x' + Δx) = y
    # A x' = y'
    # y' = y - A*Δx
    # x = A \ y' + Δx
    fill!(buf, zero(eltype(buf)))
    n = min(length(r)÷bz, length(Δp))
    pressure_indices = 1:bz:(1 + (n - 1)*bz)
    pressure_increment = view(Δp, 1:n)
    buffer_pressure = view(buf, pressure_indices)
    solution_pressure = view(x, pressure_indices)
    @tic "p increment" begin
        buffer_pressure .= pressure_increment
        solution_pressure .+= pressure_increment
    end
    correct_residual!(r, A_ps, buf)
end

@inline function set_dp!(x, bz, Δp, i)
    @inbounds x[(i-1)*bz + 1] = Δp[i]
    @inbounds for j = 2:bz
        x[(i-1)*bz + j] = 0.0
    end
end

function increment_pressure!(x, Δp, bz, p_buf = missing)
    n = min(length(x)÷bz, length(Δp))
    pressure = view(x, 1:bz:(1 + (n - 1)*bz))
    pressure .+= view(Δp, 1:n)
    return x
end


function should_update_pressure_subsystem(cpr, rec)
    interval = cpr.update_interval
    if isnothing(cpr.A_p)
        update = true
    elseif interval == :once
        update = false
    else
        it = Jutul.subiteration(rec)
        outer_step = Jutul.step(rec)
        ministep = Jutul.substep(rec)
        if interval == :iteration
            crit = true
            n = it
        elseif interval == :ministep
            n = ministep
            crit = it == 1
        elseif interval == :step
            n = outer_step
            crit = it == 1
        else
            error("Bad parameter update_frequency: $interval")
        end
        uf = cpr.update_frequency
        update = crit && (uf == 1 || (n % uf) == 1)
    end
    return update
end

function should_update_cpr(cpr, rec, type = :amg)
    if type == :partial
        interval, update_frequency = cpr.update_interval_partial, cpr.update_frequency_partial
        if !cpr.partial_update
            return false
        end
    else
        @assert type == :amg
        interval, update_frequency = cpr.update_interval, cpr.update_frequency
    end
    if isnothing(cpr.storage)
        update = true
    elseif interval == :once
        update = false
    else
        it = Jutul.subiteration(rec) + 1
        outer_step = Jutul.step(rec)
        ministep = Jutul.substep(rec)
        if interval == :iteration
            crit = true
            n = it
        elseif interval == :ministep
            n = ministep
            crit = it == 1
        elseif interval == :step
            n = outer_step
            crit = it == 1 && ministep == 1
        else
            error("Bad parameter update_frequency: $interval")
        end
        uf = update_frequency
        update = crit && (uf == 1 || (n % uf) == 1)
    end
    return update
end

function cpr_construct_well_reservoir_map(model::SimulationModel, lsys, bz)
    return nothing
end

function cprw_pattern_indices(matrix)
    rows, columns, = findnz(matrix)
    return rows, columns
end

function cprw_pattern_indices(matrix::Jutul.StaticSparsityMatrixCSR)
    rowptr = Array(matrix.rowptr)
    colval = Array(Jutul.colvals(matrix))
    rows = Vector{eltype(colval)}(undef, length(colval))
    for row in axes(matrix, 1)
        for position in rowptr[row]:(rowptr[row + 1] - 1)
            rows[position] = row
        end
    end
    return rows, colval
end

function cpr_host_matrix(matrix::Jutul.StaticSparsityMatrixCSR)
    if matrix.rowptr isa Vector
        return matrix
    end
    return Jutul.StaticSparsityMatrixCSR(
        Array(nonzeros(matrix)), Array(Jutul.colvals(matrix)),
        Array(matrix.rowptr), size(matrix, 1), size(matrix, 2), nothing;
        nthreads = 1, minbatch = matrix.minbatch,
        thread_type = :serial)
end

cpr_host_matrix(matrix) = matrix

function refresh_cprw_source_values!(map)
    sources = map.source_values
    if !isnothing(sources)
        for (destination, blocks) in zip((map.nzval_12, map.nzval_21, map.nzval_22), sources)
            offset = 0
            for block in blocks
                copyto!(destination, offset + 1, block, 1, length(block))
                offset += length(block)
            end
        end
    end
    return map
end

function cpr_construct_well_reservoir_map(model::MultiModel, lsys, bz)
    num_wells = 0
    well_keys = Symbol[]
    well_ranges = UnitRange{Int}[]
    well_nodes = Vector{UnitRange{Int}}[]
    for (k, m) in pairs(model.models)
        if !model_or_domain_is_well(m)
            continue
        end
        w = physical_representation(m.domain)
        nodes = if w isa SimpleWell
            [node:node for node in 1:number_of_cells(w)]
        elseif w isa MultiSegmentWell
            # A uniform pressure shift over a physical well provides its BHP
            # coarse mode; segment fluxes and transport stay in the Schur solve.
            if isnothing(w.multiwell)
                [1:number_of_cells(w)]
            else
                [w.multiwell.nodes.pos[i]:(w.multiwell.nodes.pos[i + 1] - 1)
                    for i in eachindex(w.multiwell.names)]
            end
        else
            continue
        end
        n = length(nodes)
        push!(well_ranges, (num_wells + 1):(num_wells + n))
        push!(well_nodes, nodes)
        num_wells += n
        push!(well_keys, k)
    end
    isempty(well_keys) && return nothing

    # Split wells can belong to different Schur groups. Keep direct references
    # for the usual single group; multiple groups share a refreshed flat buffer
    # so the numerical reduction and backend kernels use the same mappings.
    groups = unique([Jutul.group_index(model, key) for key in well_keys])
    sources = map(((1, 2), (2, 1), (2, 2))) do (row, col)
        [nonzeros(lsys[row == 1 ? 1 : group, col == 1 ? 1 : group].jac) for group in groups]
    end
    offsets = map(blocks -> cumsum([0; length.(blocks)]), sources)
    values = map(blocks -> length(blocks) == 1 ? only(blocks) : reduce(vcat, blocks), sources)

    jacobian_layout(::Jutul.LinearizedBlock{R}) where R = R()
    jacobian_layout(system::LinearizedSystem) = system.matrix_layout
    function well_positions(wcells, wperf, system, row_index, column_index,
            ::Val{bz}, nz_offset) where bz
        jacobian = cpr_host_matrix(system.jac)
        adjoint = Jutul.represented_as_adjoint(jacobian_layout(system))
        return map(eachindex(wcells)) do i
            out = @MMatrix zeros(Int64, bz, bz)
            for der in 1:bz, eq in 1:bz
                row = row_index(wcells[i], eq)
                column = column_index(wperf[i], der)
                # Multisegment conservation blocks can have structural zeros:
                # segment fluxes are independent unknowns on faces.
                pos = find_sparse_position(jacobian, row, column, adjoint)
                out[eq, der] = iszero(pos) ? 0 : pos + nz_offset
            end
            return SMatrix{bz, bz, Int64, bz*bz}(out)
        end
    end

    # For each well find the list of cells and corresponding sparsity + indices
    # since these are not blocks.
    J_rr = lsys[1, 1].jac
    I, J = cprw_pattern_indices(J_rr)
    ncell = size(J_rr, 1)

    bz = bz + model_is_thermal(reservoir_model(model))
    reservoir_index(cell, component) = bz*(cell - 1) + component
    well_cells = Vector{Int64}[]
    map_T = Vector{Vector{SMatrix{bz, bz, Int64, bz*bz}}}
    map_12 = map_T()
    map_21 = map_T()
    map_22 = map_T()

    for (group, k) in enumerate(well_keys)
        wmodel = model.models[k]
        w = physical_representation(wmodel.domain)
        nwell_cells = number_of_cells(w)
        linear_group = Jutul.group_index(model, k)
        source_index = findfirst(isequal(linear_group), groups)
        model_offset = Jutul.group_linearized_system_offset(model, k)
        well_layout = matrix_layout(wmodel.context)
        variable_offset = model_offset
        for entity in Jutul.get_primary_variable_ordered_entities(wmodel)
            entity == Cells() && break
            variable_offset += Jutul.number_of_degrees_of_freedom(wmodel, entity)
        end
        equation_offset = Jutul.group_linearized_system_offset(model, k, Jutul.number_of_equations) +
            Jutul.get_equation_offset(wmodel, :mass_conservation)
        well_column(cell, component) = variable_offset +
            Jutul.alignment_linear_index(cell, component, nwell_cells, bz, well_layout)
        well_row(cell, component) = equation_offset +
            Jutul.alignment_linear_index(cell, component, nwell_cells, bz, well_layout)
        # This map is constructed once on the host. The completed index arrays
        # are moved to the execution backend by prepare_cprw_backend_maps!.
        perforations = w.perforations
        connections = unique(collect(zip(Array(perforations.self), Array(perforations.reservoir))))
        for (coarse_node, i) in enumerate(well_ranges[group])
            nodes = well_nodes[group][coarse_node]
            perfs = filter(pair -> first(pair) in nodes, connections)
            wn = first.(perfs)
            wc = last.(perfs)
            push!(well_cells, wc)
            # Reservoir differentiated with respect to wells
            push!(map_12, well_positions(wc, wn,
                lsys[1, linear_group], reservoir_index, well_column, Val(bz), offsets[1][source_index]))
            # Wells differentiated with respect to reservoir
            push!(map_21, well_positions(wn, wc,
                lsys[linear_group, 1], well_row, reservoir_index, Val(bz), offsets[2][source_index]))
            # Wells differentiated with respect to wells
            diagonal_pairs = [(node, node) for node in nodes]
            if w isa MultiSegmentWell
                for (a, b) in eachcol(Array(w.neighborship))
                    if a in nodes && b in nodes
                        push!(diagonal_pairs, (a, b), (b, a))
                    end
                end
            end
            unique!(diagonal_pairs)
            push!(map_22, well_positions(first.(diagonal_pairs), last.(diagonal_pairs),
                lsys[linear_group, linear_group], well_row, well_column, Val(bz), offsets[3][source_index]))
            # Extend sparsity to include well cells
            for c in wc
                push!(I, c)
                push!(J, i + ncell)

                push!(J, c)
                push!(I, i + ncell)
            end
            push!(I, i + ncell)
            push!(J, i + ncell)
        end
    end

    nzmap_reservoir = Int[]
    nzmap_12 = Int[]
    nzmap_21 = Int[]
    nzmap_well = Int[]

    well_reservoir_map = (
        wells = well_keys,
        well_ranges = well_ranges,
        I = I,
        J = J,
        np = ncell + num_wells,
        well_cells = well_cells,
        map_12 = map_12,
        map_21 = map_21,
        map_22 = map_22,
        nzmap_11 = nzmap_reservoir,
        nzmap_12 = nzmap_12,
        nzmap_21 = nzmap_21,
        nzmap_22 = nzmap_well,
        nzval_12 = values[1],
        nzval_21 = values[2],
        nzval_22 = values[3],
        source_values = length(groups) == 1 ? nothing : sources,
        backend_data = Ref{Any}(nothing)
    )
    return well_reservoir_map
end
