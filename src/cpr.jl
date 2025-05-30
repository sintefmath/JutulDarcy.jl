
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
    solution = zeros(T, ncell*bz)
    residual = zeros(T, ncell*bz)
    w_p = zeros(T, ncomp, np)
    w_rhs = zeros(ncomp)
    w_rhs[1] = 1
    w_rhs = SVector{ncomp, T}(w_rhs)
    if p_buffer
        p_buf = zeros(np)
    else
        p_buf = missing
    end
    return CPRStorage(A_p, r_p, p, solution, residual, lin_op, w_p, w_rhs, np, bz, ncomp, objectid(full_jac), p_buf, well_reservoir_map, T)
end

function CPRStorage(np::Int, bz::Int, lin_op, psys::Tuple, solution, residual, T = Float64, id = zero(UInt64); ncomp = bz)
    A_p, r_p, p = psys
    w_p = zeros(T, ncomp, np)
    w_rhs = zeros(ncomp)
    w_rhs[1] = 1
    w_rhs = SVector{bz, T}(w_rhs)
    return CPRStorage(A_p, r_p, p, solution, residual, lin_op, w_p, w_rhs, np, bz, ncomp, id, zeros(T, np), nothing, T)
end


"""
    CPRPreconditioner(p = default_psolve(), s = ILUZeroPreconditioner(); strategy = :quasi_impes, weight_scaling = :unit, update_frequency = 1, update_interval = :iteration, partial_update = true)

Construct a constrained pressure residual (CPR) preconditioner.

By default, this is a AMG-BILU(0) version (algebraic multigrid for pressure,
block-ILU(0) for the global system).
"""
mutable struct CPRPreconditioner{P, S} <: JutulPreconditioner
    storage::Union{CPRStorage, Nothing}
    pressure_precond::P
    system_precond::S
    strategy::Symbol
    variant::Symbol
    weight_scaling::Symbol
    update_frequency::Int             # Update frequency for AMG hierarchy (and pressure part if partial_update = false)
    update_interval::Symbol           # iteration, ministep, step, ...
    update_frequency_partial::Int     # Update frequency for pressure system
    update_interval_partial::Symbol   # iteration, ministep, step, ...
    partial_update::Bool              # Perform partial update of AMG and update pressure system
    full_system_correction::Bool
    p_rtol::Union{Float64, Nothing}
    npre::Int
    npost::Int
    mode::Symbol
    psolver
end

function CPRPreconditioner(p = default_psolve(), s = ILUZeroPreconditioner();
        strategy = :true_impes,
        variant = :cpr,
        weight_scaling = :unit,
        npre = 0,
        npost = 1,
        update_frequency = 1,
        update_interval = :iteration,
        update_frequency_partial = 1,
        update_interval_partial = :iteration,
        p_rtol = nothing,
        full_system_correction = true,
        partial_update = update_interval == :once,
        mode = :forward
    )
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
        full_system_correction,
        p_rtol,
        npre,
        npost,
        mode,
        nothing
    )
end

function default_psolve(; max_levels = 10, max_coarse = 10, amgcl_type = :amg, type = default_amg_symbol(), kwarg...)
    if type == :hypre
        amg = BoomerAMGPreconditioner(; kwarg...)
    elseif type == :amgx
        amg = AMGXPreconditioner(; kwarg...)
    elseif type == :amgcl
        if length(kwarg) == 0
            # Some reasonable defaults for reservoir system
            agg = (
                coarsening = (
                    type = "aggregation",
                    over_interp = 1.0,
                    aggr = (
                        eps_strong = 0.1,
                    )
                ),
                npre = 3,
                npost = 3,
                ncycle = 1,
                coarse_enough = 1000,
                pre_cycles = 1,
                relax = (
                    type = "spai0",
                ),
            )
            if amgcl_type == :amg
                # Direct AMG as preconditioner
                kwarg = agg
            elseif amgcl_type == :amg_solver
                # Nexted Krylov solve - should use FGMRES on outside.
                kwarg = (
                    solver = (
                        type = :fgmres,
                        tol = 1e-2,
                        verbose = false
                        ),
                    precond = agg,
                )
            end
        end
        amg = Jutul.AMGCLPreconditioner(amgcl_type; kwarg...)
    else
        gs_its = 1
        cyc = AlgebraicMultigrid.V()
        gs = GaussSeidel(ForwardSweep(), gs_its)
        amg = AMGPreconditioner(type; max_levels = max_levels, max_coarse = max_coarse, presmoother = gs, postsmoother = gs, cycle = cyc, kwarg...)
    end
end

function update_preconditioner!(cpr::CPRPreconditioner, lsys::Jutul.JutulLinearSystem, ctx_outer, model, storage, recorder, executor; update_system_precond = true, T = Float64)
    rmodel = reservoir_model(model, type = :flow)
    ctx = rmodel.context
    update_p = update_cpr_internals!(cpr, lsys, model, storage, recorder, executor, T)
    if update_system_precond
        @tic "s-precond" update_preconditioner!(cpr.system_precond, lsys, ctx, model, storage, recorder, executor)
    end
    if update_p
        @tic "p-precond" update_preconditioner!(cpr.pressure_precond, cpr.storage.A_p, cpr.storage.r_p, ctx, executor)
    elseif should_update_cpr(cpr, recorder, :partial)
        @tic "p-precond (partial)" partial_update_preconditioner!(cpr.pressure_precond, cpr.storage.A_p, cpr.storage.r_p, ctx, executor)
    end
end

function initialize_cpr_storage!(cpr, model, lsys, s, bz, T)
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
    end
end

function pressure_matrix_from_global_jacobian(J::SparseMatrixCSC, T, lsys, well_reservoir_map::Nothing)
    nzval = zeros(T, nnz(J))
    n = size(J, 2)
    return SparseMatrixCSC(n, n, J.colptr, J.rowval, nzval)
end

function pressure_matrix_from_global_jacobian(J::Jutul.StaticSparsityMatrixCSR, T, lsys, well_reservoir_map::Nothing)
    nzval = zeros(T, nnz(J))
    n = size(J, 2)
    # Assume symmetry in sparse pattern, but not values.
    return Jutul.StaticSparsityMatrixCSR(n, n, J.At.colptr, Jutul.colvals(J), nzval, nthreads = J.nthreads, minbatch = J.minbatch)
end

function pressure_matrix_from_global_jacobian(sys_jac::Jutul.StaticSparsityMatrixCSR, T, lsys, well_reservoir_map)
    n = well_reservoir_map.np
    I = well_reservoir_map.I
    J = well_reservoir_map.J
    @assert length(I) == length(J)
    V = zeros(T, length(I))
    n::Integer
    A = sparse(J, I, V, n, n)
    A_p = Jutul.StaticSparsityMatrixCSR(A;
        nthreads = sys_jac.nthreads,
        minbatch = sys_jac.minbatch
    )
    # Align reservoir variables
    nzmap_reservoir = well_reservoir_map.nzmap_11
    resize!(nzmap_reservoir, nnz(sys_jac))
    ix_in_pnzval = 1
    ncell = size(sys_jac, 1)
    nwell = size(A_p, 1) - ncell
    for row in 1:ncell
        for j in nzrange(sys_jac, row)
            col = Jutul.StaticCSR.colvals(sys_jac)[j]
            nzmap_reservoir[j] = find_sparse_position(A_p, row, col)
        end
    end
    # Well to well
    nzmap_well = well_reservoir_map.nzmap_22
    for w in 1:nwell
        ix = find_sparse_position(A_p, ncell + w, ncell + w)
        push!(nzmap_well, ix)
    end

    # Align well to reservoir variables
    nzmap_12 = well_reservoir_map.nzmap_12
    nzmap_21 = well_reservoir_map.nzmap_21
    J_12 = lsys[1, 2].jac
    J_21 = lsys[2, 1].jac
    for w in 1:nwell
        for (i, cell) in enumerate(well_reservoir_map.well_cells[w])
            ix_12 = find_sparse_position(A_p, cell, ncell + w)
            push!(nzmap_12, ix_12)
            ix_21 = find_sparse_position(A_p, ncell + w, cell)
            push!(nzmap_21, ix_21)
        end
    end
    return A_p
end

function create_pressure_system(p_prec, J, lsys, n, T, well_reservoir_map)
    J_p = pressure_matrix_from_global_jacobian(J, T, lsys, well_reservoir_map)
    return (J_p, zeros(T, n), zeros(T, n))
end

function update_cpr_internals!(cpr::CPRPreconditioner, lsys, model, storage, recorder, executor, T)
    do_p_update = should_update_cpr(cpr, recorder, :amg)
    s = reservoir_storage(model, storage)
    A = reservoir_jacobian(lsys)
    rmodel = reservoir_model(model, type = :flow)
    bz = number_of_components(rmodel.system)
    initialize_cpr_storage!(cpr, model, lsys, s, bz, T)
    ps = rmodel.primary_variables[:Pressure].scale
    if do_p_update || cpr.partial_update
        rmodel = reservoir_model(model)
        ctx = rmodel.context
        @tic "weights" w_p = update_weights!(cpr, cpr.storage, rmodel, s, storage, A, ps)
        A_p = cpr.storage.A_p
        w_p = cpr.storage.w_p
        rw_map = cpr.storage.well_reservoir_map
        @tic "pressure system" update_pressure_system!(A_p, cpr.pressure_precond, A, w_p, ctx, executor, rw_map)
    end
    return do_p_update
end

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
            Ji = nz_12[map_12_i[j]]
            nz_p[nzmap_12[ix]] = reduce_to_pressure(Ji, w_p, cell, ncomp, is_adjoint)
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
            Ji = nz_21[map_21_i[j]]
            nz_p[nzmap_21[ix]] = reduce_to_pressure(Ji, w_p, well, ncomp, is_adjoint)
            ix += 1
        end
    end
    # Wells differentiated wrt wells
    nz_22 = well_reservoir_map.nzval_22
    ix = 1
    @inbounds for i in eachindex(map_22, nzmap_22)
        map_22_i = map_22[i]
        @inbounds for j in eachindex(map_22_i)
            well = i + ncells
            Ji = nz_22[map_22_i[j]]
            nz_p[nzmap_22[ix]] = reduce_to_pressure(Ji, w_p, well, ncomp, is_adjoint)
            ix += 1
        end
    end
    return nz_p
end

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
    @. x = 0.0
    # presmooth
    @tic "cpr smoother" apply_cpr_smoother!(x, r, buf, smoother, A_ps, cpr.npre)
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
    @tic "residual correction" mul!(r, A, x, -1.0, true)
end

function apply_cpr_pressure_stage!(cpr::CPRPreconditioner, cpr_s::CPRStorage, r, arg...)
    r_p, w_p, bz, Δp = cpr_s.r_p, cpr_s.w_p, cpr_s.block_size, cpr_s.p
    ncomp = cpr_s.number_of_components
    @tic "p rhs" update_p_rhs!(r_p, r, ncomp, bz, w_p, cpr.pressure_precond, cpr.mode)
    # Apply preconditioner to pressure part
    @tic "p apply" begin
        p_rtol = cpr.p_rtol
        p_precond = cpr.pressure_precond
        cpr_p_apply!(Δp, cpr, p_precond, r_p, p_rtol)
    end
end

function cpr_p_apply!(Δp, cpr, p_precond, r_p, p_rtol)
    apply!(Δp, p_precond, r_p)
    if !isnothing(p_rtol)
        A_p = cpr.A_p
        if isnothing(cpr.psolver)
            cpr.psolver = FgmresSolver(A_p, r_p)
        end
        psolve = cpr.psolver
        warm_start!(psolve, Δp)
        M = Jutul.PrecondWrapper(linear_operator(p_precond))
        fgmres!(psolve, A_p, r_p, M = M, rtol = p_rtol, atol = 1e-12, itmax = 20)
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

function update_weights!(cpr, cpr_storage::CPRStorage, model, res_storage, storage, J, ps)
    np = cpr_storage.np
    nc = number_of_cells(reservoir_model(model).domain)
    n = min(np, nc)
    bz = cpr_storage.block_size
    w = cpr_storage.w_p
    r = cpr_storage.w_rhs
    ncomp = size(w, 1)
    scaling = cpr.weight_scaling
    strategy = cpr.strategy
    if strategy == :true_impes
        eq_s = res_storage.equations[:mass_conservation]
        if eq_s isa ConservationLawTPFAStorage
            acc = eq_s.accumulation.entries
        else
            acc = res_storage.state.TotalMasses
            # This term isn't scaled by dt, so use simple weights instead
            ps = 1.0
        end
        true_impes!(w, acc, r, n, ncomp, ps, scaling)
    elseif strategy == :analytical
        rstate = res_storage.state
        cpr_weights_no_partials!(w, model, rstate, r, n, ncomp, scaling)
    elseif strategy == :quasi_impes
        quasi_impes!(w, J, r, n, ncomp, scaling)
    elseif strategy == :none
        # Do nothing. Already set to one.
    else
        error("Unsupported strategy $(strategy)")
    end
    if np > nc
        wr_map = cpr_storage.well_reservoir_map
        @assert !isnothing(wr_map)
        for i in 1:(np-nc)
            wno = nc+i
            w_i = view(w, :, wno)
            well = wr_map.wells[i]
            wstate = storage[well].state
            if strategy == :analytical
                cpr_weights_no_partials!(w_i, model, wstate, nothing, 1, ncomp, scaling)
            else
                acc_i = wstate.TotalMasses
                true_impes!(w_i, acc_i, r, 1, ncomp, ps, scaling)
            end
        end
    end
    return w
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

function quasi_impes!(w, J, r, n, bz, scaling)
    r_p = SVector{bz}(r)
    @batch for cell = 1:n
        J_b = J[cell, cell]'
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
    @. buf = 0
    n = length(r)÷bz
    @tic "p increment" @inbounds for i in 1:n
        ix = (i-1)*bz + 1
        p_i = Δp[i]
        buf[ix] = p_i
        x[ix] += p_i
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
    @inbounds for i in eachindex(Δp)
        x[(i-1)*bz + 1] += Δp[i]
    end
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

function cpr_construct_well_reservoir_map(model::MultiModel, lsys, bz)
    num_simple_wells = 0
    simple_well_keys = Symbol[]
    for (k, m) in pairs(model.models)
        if k == :Reservoir
            continue
        end
        stop = true
        if model_or_domain_is_well(m)
            if physical_representation(m.domain) isa SimpleWell
                num_simple_wells += 1
                push!(simple_well_keys, k)
                stop = false
            end
        end
        if stop
            break
        end
    end

    function weq_positions(J, cell_i, cell_j, ::Val{bz}, layout) where bz
        out = @MMatrix zeros(Int64, bz, bz)
        for eq in 1:bz
            for der in 1:bz
                out[eq, der] = find_sparse_position(J, cell_i+eq-1, cell_j+der-1, layout)
            end
        end
        return SMatrix{bz, bz, Int64, bz*bz}(out)
    end

    function well_positions(wcells, wperf, lsys_block::Jutul.LinearizedBlock{R}, bz) where R
        vbz = Val(bz)
        layout = R()
        return map(i -> weq_positions(lsys_block.jac, bz*(wcells[i]-1)+1, bz*(wperf[i]-1)+1, vbz, layout), eachindex(wcells))
    end

    function well_positions(wcells, wperf, lsys::LinearizedSystem, bz)
        vbz = Val(bz)
        layout = lsys.matrix_layout
        return map(i -> weq_positions(lsys.jac, bz*(wcells[i]-1)+1, bz*(wperf[i]-1)+1, vbz, layout), eachindex(wcells))
    end


    # For each well find the list of cells and corresponding sparsity + indices
    # since these are not blocks.
    J_rr = lsys[1, 1].jac
    if J_rr isa Jutul.StaticSparsityMatrixCSR
        # TODO: Add this method to CSR
        J, I, = findnz(J_rr.At)
    else
        I, J, = findnz(J_rr)
    end
    ncell = size(J_rr, 1)

    well_cells = Vector{Int64}[]
    map_T = Vector{Vector{SMatrix{bz, bz, Int64, bz*bz}}}
    map_12 = map_T()
    map_21 = map_T()
    map_22 = map_T()
    map_well = map_T()

    for (i, k) in enumerate(simple_well_keys)
        wmodel = model.models[k]
        # layout = matrix_layout(wmodel.context)
        w = physical_representation(wmodel.domain)
        wc = w.perforations.reservoir
        push!(well_cells, wc)
        vbz = Val(bz)

        wno = fill(i, length(wc))
        # Reservoir differentiated with respect to wells
        pos_12 = well_positions(wc, wno, lsys[1, 2], bz)
        push!(map_12, pos_12)
        # Wells differentiated with respect to reservoir
        pos_21 = well_positions(wno, wc, lsys[2, 1], bz)
        push!(map_21, pos_21)
        # Wells differentiated with respect to wells
        pos_22 = well_positions([i], [i], lsys[2, 2], bz)
        push!(map_22, pos_22)
        # Extend sparsity to include well cells
        for c in wc
            push!(I, c)
            push!(J, i + ncell)

            push!(J, c)
            push!(I, i + ncell)

            push!(I, i + ncell)
            push!(J, i + ncell)
        end
    end

    nzmap_reservoir = Int[]
    nzmap_12 = Int[]
    nzmap_21 = Int[]
    nzmap_well = Int[]

    well_reservoir_map = (
        wells = simple_well_keys,
        I = I,
        J = J,
        np = ncell + num_simple_wells,
        well_cells = well_cells,
        map_12 = map_12,
        map_21 = map_21,
        map_22 = map_22,
        nzmap_11 = nzmap_reservoir,
        nzmap_12 = nzmap_12,
        nzmap_21 = nzmap_21,
        nzmap_22 = nzmap_well,
        nzval_12 = nonzeros(lsys[1, 2].jac),
        nzval_21 = nonzeros(lsys[2, 1].jac),
        nzval_22 = nonzeros(lsys[2, 2].jac)
    )
    return well_reservoir_map
end
