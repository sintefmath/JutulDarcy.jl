export CPRPreconditioner
"""
Constrained pressure residual
"""
mutable struct CPRPreconditioner <: JutulPreconditioner
    A_p  # pressure system
    r_p  # pressure residual
    p    # last pressure approximation
    buf  # buffer of size equal to full system rhs
    A_ps # full system
    w_p  # pressure weights
    pressure_precond
    system_precond
    strategy
    weight_scaling
    block_size
    update_frequency::Integer # Update frequency for AMG hierarchy (and pressure part if partial_update = false)
    update_interval::Symbol   # iteration, ministep, step, ...
    partial_update            # Perform partial update of AMG and update pressure system
    function CPRPreconditioner(p = default_psolve(), s = ILUZeroPreconditioner(); strategy = :quasi_impes, weight_scaling = :unit, update_frequency = 1, update_interval = :iteration, partial_update = true)
        new(nothing, nothing, nothing, nothing, nothing, nothing, p, s, strategy, weight_scaling, nothing, update_frequency, update_interval, partial_update)
    end
end

function default_psolve(; max_levels = 10, max_coarse = 10, type = :smoothed_aggregation)
    gs_its = 1
    cyc = AlgebraicMultigrid.V()
    if type == :smoothed_aggregation
        m = smoothed_aggregation
    else
        m = ruge_stuben
    end
    gs = GaussSeidel(ForwardSweep(), gs_its)
    return AMGPreconditioner(m, max_levels = max_levels, max_coarse = max_coarse, presmoother = gs, postsmoother = gs, cycle = cyc)
end

function update!(cpr::CPRPreconditioner, lsys, model, arg...)
    rmodel = reservoir_model(model)
    ctx = rmodel.context
    update_p = update_cpr_internals!(cpr, lsys, model, arg...)
    @timeit "s-precond" update!(cpr.system_precond, lsys, model, arg...)
    if update_p
        @timeit "p-precond" update!(cpr.pressure_precond, cpr.A_p, cpr.r_p, ctx)
    elseif cpr.partial_update
        @timeit "p-precond (partial)" partial_update!(cpr.pressure_precond, cpr.A_p, cpr.r_p, ctx)
    end
end

function initialize_storage!(cpr, J, s)
    if isnothing(cpr.A_p)
        m, n = size(J)
        cpr.block_size = bz = size(eltype(J), 1)
        # @assert n == m == length(s.state.Pressure) "Expected Jacobian dimensions ($m by $n) to both equal number of pressures $(length(s.state.Pressure))"
        cpr.A_p = create_pressure_matrix(J)
        cpr.r_p = zeros(n)
        cpr.buf = zeros(n*bz)
        cpr.p = zeros(n)
        cpr.w_p = zeros(bz, n)
    end
end

function create_pressure_matrix(J)
    nzval = zeros(nnz(J))
    n = size(J, 2)
    SparseMatrixCSC(n, n, J.colptr, J.rowval, nzval)
end


function create_pressure_matrix(J::Jutul.StaticSparsityMatrixCSR)
    nzval = zeros(nnz(J))
    n = size(J, 2)
    # Assume symmetry in sparse pattern, but not values.
    Jutul.StaticSparsityMatrixCSR(n, n, J.At.colptr, Jutul.colvals(J), nzval, nthreads = J.nthreads, minbatch = J.minbatch)
end

function update_cpr_internals!(cpr::CPRPreconditioner, lsys, model, storage, recorder)
    do_p_update = should_update_pressure_subsystem(cpr, recorder)
    s = reservoir_storage(model, storage)
    A = reservoir_jacobian(lsys)
    rmodel = reservoir_model(model)
    cpr.A_ps = linear_operator(lsys)
    initialize_storage!(cpr, A, s)
    ps = rmodel.primary_variables[:Pressure].scale
    if do_p_update || cpr.partial_update
        @timeit "weights" w_p = update_weights!(cpr, rmodel, s, A, ps)
        @timeit "pressure system" update_pressure_system!(cpr.A_p, A, w_p, cpr.block_size, model.context)
    end
    return do_p_update
end

function update_pressure_system!(A_p, A, w_p, bz, ctx)
    cp = A_p.colptr
    nz = A_p.nzval
    nz_s = A.nzval
    rv = A_p.rowval
    @assert size(nz) == size(nz_s)
    n = A.n
    # Update the pressure system with the same pattern in-place
    tb = minbatch(ctx)
    @batch minbatch=tb for i in 1:n
        @inbounds for j in cp[i]:cp[i+1]-1
            row = rv[j]
            Ji = nz_s[j]
            tmp = 0.0
            @inbounds for b = 1:bz
                tmp += Ji[b, 1]*w_p[b, row]
            end
            nz[j] = tmp
        end
    end
end

function update_pressure_system!(A_p::Jutul.StaticSparsityMatrixCSR, A::Jutul.StaticSparsityMatrixCSR, w_p, bz, ctx)
    nz = nonzeros(A_p)
    nz_s = nonzeros(A)
    cols = Jutul.colvals(A)
    @assert size(nz) == size(nz_s)
    n = size(A_p, 1)
    # Update the pressure system with the same pattern in-place
    tb = minbatch(ctx)
    @batch minbatch=tb for i in 1:n
        @inbounds for j in nzrange(A, i)
            col = cols[j]
            Ji = nz_s[j]
            tmp = 0.0
            @inbounds for b = 1:bz
                tmp += Ji[b, 1]*w_p[b, col]
            end
            nz[j] = tmp
        end
    end
end

function operator_nrows(cpr::CPRPreconditioner)
    return length(cpr.r_p)*cpr.block_size
end

function apply!(x, cpr::CPRPreconditioner, r, arg...)
    r_p, w_p, bz, Δp = cpr.r_p, cpr.w_p, cpr.block_size, cpr.p
    if false
        y = copy(r)
        # y = r
        # Construct right hand side by the weights
        norm0 = norm(r)
        do_cpr = true
        update_p_rhs!(r_p, y, bz, w_p)
        println("**************************************************************")
        # Apply preconditioner to pressure part
        @info "Before pressure correction" norm(y) norm(r_p)
        if do_cpr
            apply!(Δp, cpr.pressure_precond, r_p)
            correct_residual_for_dp!(y, x, Δp, bz, cpr.buf, cpr.A_ps)
            norm_after = norm(y)
            @info "After pressure correction" norm(y) norm(cpr.A_p*Δp - r_p) norm_after/norm0
        end
        apply!(x, cpr.system_precond, y)
        if do_cpr
            @info "After second stage" norm(cpr.A_ps*x - y)
            increment_pressure!(x, Δp, bz)
        end
        @info "Final" norm(cpr.A_ps*x - r) norm(cpr.A_ps*x - r)/norm0
    else
        y = r
        # Construct right hand side by the weights
        @timeit "p rhs" update_p_rhs!(r_p, y, bz, w_p)
        # Apply preconditioner to pressure part
        @timeit "p apply" apply!(Δp, cpr.pressure_precond, r_p)
        @timeit "r update" correct_residual_for_dp!(y, x, Δp, bz, cpr.buf, cpr.A_ps)
        @timeit "s apply" apply!(x, cpr.system_precond, y)
        @timeit "Δp" increment_pressure!(x, Δp, bz)
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

function update_weights!(cpr, model, res_storage, J, ps)
    n = size(cpr.A_p, 1)
    bz = cpr.block_size
    if isnothing(cpr.w_p)
        cpr.w_p = ones(bz, n)
    end
    w = cpr.w_p
    r = zeros(bz)
    r[1] = 1.0
    scaling = cpr.weight_scaling
    if cpr.strategy == :true_impes
        eq = res_storage.equations[:mass_conservation]
        acc = eq.accumulation.entries
        true_impes!(w, acc, r, n, bz, ps, scaling)
    elseif cpr.strategy == :analytical
        rstate = res_storage.state
        cpr_weights_no_partials!(w, model, rstate, r, n, bz, scaling)
    elseif cpr.strategy == :quasi_impes
        quasi_impes!(w, J, r, n, bz, scaling)
    elseif cpr.strategy == :none
        # Do nothing. Already set to one.
    else
        error("Unsupported strategy $(cpr.strategy)")
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
        true_impes_gen!(w, acc, r, n, bz, arg...)
    end
end

function true_impes_2!(w, acc, r, n, bz, p_scale, scaling)
    r_p = SVector{2}(r)
    for cell in 1:n
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

function true_impes_gen!(w, acc, r, n, bz, p_scale, scaling)
    r_p = SVector{bz}(r)
    A = MMatrix{bz, bz, eltype(r)}(zeros(bz, bz))
    for cell in 1:n
        @inbounds for i = 1:bz
            v = acc[i, cell]
            @inbounds A[1, i] = v.partials[1]*p_scale
            @inbounds for j = 2:bz
                A[j, i] = v.partials[j]
            end
        end
        invert_w!(w, A, r_p, cell, bz, scaling)
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
        s = 1.0/sum(tmp)
    else
        s = 1.0
    end
    @inbounds for i = 1:bz
        w[i, cell] = tmp[i]*s
    end
end

function update_p_rhs!(r_p, y, bz, w_p)
    if false
        @batch minbatch = 1000 for i in eachindex(r_p)
            v = 0.0
            @inbounds for b = 1:bz
                v += y[(i-1)*bz + b]*w_p[b, i]
            end
            @inbounds r_p[i] = v
        end
    end
    n = length(y) ÷ bz
    yv = reshape(y, bz, n)
    @tullio r_p[i] = yv[b, i]*w_p[b, i]
end

function correct_residual_for_dp!(y, x, Δp, bz, buf, A)
    # x = x' + Δx
    # A (x' + Δx) = y
    # A x' = y'
    # y' = y - A*Δx
    # x = A \ y' + Δx
    @batch minbatch = 1000 for i in eachindex(Δp)
        set_dp!(x, bz, Δp, i)
    end
    if false
        mul!(buf, A, x)
        @batch minbatch = 1000 for i in eachindex(y)
            @inbounds y[i] -= buf[i]
        end
    else
        mul!(y, A, x, -1.0, 1.0)
    end
    # @. y -= buf
end

@inline function set_dp!(x, bz, Δp, i)
    @inbounds x[(i-1)*bz + 1] = Δp[i]
    @inbounds for j = 2:bz
        x[(i-1)*bz + j] = 0.0
    end
end

function increment_pressure!(x, Δp, bz)
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