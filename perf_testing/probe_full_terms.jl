# Additional fixed-input terms: full well-reduced operator, CPR transfers and
# the individual parts of an AMG numeric update. Reuse the measurement helper.
include("probe_terms.jl")
using SparseArrays

function probe_full_terms()
    fixture = deserialize(ARGS[1])
    cpr, lsys = fixture.preconditioner, fixture.linear_system
    rhs = copy(cpr.storage.r_ps)
    rhs ./= norm(rhs)
    x = zeros(eltype(rhs), length(rhs))
    prhs = copy(cpr.storage.r_p)
    prhs ./= norm(prhs)
    px = similar(prhs)
    terms, couplings, updates = Any[], Any[], Any[]
    push!(terms, term_row("saved_cpr_apply", () -> Jutul.apply!(x, cpr, rhs)))
    push!(terms, term_row("saved_pressure_apply", () -> Jutul.apply!(px, cpr.pressure_precond, prhs)))
    matrix = JutulDarcy.reservoir_jacobian(lsys)
    serial_matrix = Jutul.StaticSparsityMatrixCSR(matrix.nzval, matrix.colval, matrix.rowptr,
        size(matrix, 1), size(matrix, 2), matrix.backend;
        nthreads = 1, minbatch = matrix.minbatch, thread_type = :serial)
    native = ILUZeroPreconditioner()
    Jutul.update_preconditioner!(native, serial_matrix, rhs, nothing, nothing)
    matched_level = Jutul.KASmootherPreconditioner(:ilu0)
    Jutul.update_preconditioner!(matched_level, serial_matrix, rhs, nothing, nothing)
    reference = similar(rhs)
    Jutul.apply!(reference, matched_level, rhs)
    Jutul.apply!(x, native, rhs)
    @assert x ≈ reference
    push!(terms, term_row("native_serial_system_apply", () -> Jutul.apply!(x, native, rhs)))
    push!(terms, term_row("matched_level_system_apply", () -> Jutul.apply!(x, matched_level, rhs)))
    push!(terms, term_row("native_serial_system_update", () ->
        Jutul.update_preconditioner!(native, serial_matrix, rhs, nothing, nothing); inner = 3))
    push!(terms, term_row("full_operator", () -> mul!(x, cpr.storage.A_ps, rhs, 1.0, 0.0)))
    push!(terms, term_row("cpr_pressure_rhs", () -> JutulDarcy.update_p_rhs!(
        cpr.storage.r_p, rhs, cpr.storage.number_of_components, cpr.storage.block_size,
        cpr.storage.w_p, cpr.storage.A_p, cpr.mode, cpr.pressure_precond)))
    push!(terms, term_row("cpr_pressure_increment", () -> JutulDarcy.increment_pressure!(
        x, cpr.storage.p, cpr.storage.block_size, cpr.storage.p_buffer)))
    push!(terms, term_row("vector_copy", () -> copyto!(x, rhs); inner = 50))
    push!(terms, term_row("vector_clear", () -> fill!(x, 0.0); inner = 50))
    if lsys isa Jutul.MultiLinearizedSystem && Jutul.do_schur(lsys)
        B, C, D, E = Jutul.get_schur_blocks!(lsys, false)
        for i in eachindex(D)
            b1, b2 = lsys.schur_buffer[i + 1]
            Di, Ci, Ei = D[i], C[i], E[i]
            row = Dict("index" => i, "D_type" => string(typeof(Di)), "C_type" => string(typeof(Ci)),
                "D_size" => collect(size(Di)), "C_size" => collect(size(Ci)),
                "D_nonzeros" => nnz(Di), "C_nonzeros" => nnz(Ci), "terms" => Any[])
            push!(row["terms"], term_row("schur$(i)_D_product", () -> mul!(b2, Di, rhs)))
            push!(row["terms"], term_row("schur$(i)_E_solve", () -> ldiv!(b1, Ei, b2)))
            push!(row["terms"], term_row("schur$(i)_C_product", () -> mul!(x, Ci, b1, -1.0, 1.0)))
            push!(couplings, row)
        end
    end
    op, precond = cpr.storage.A_ps, Jutul.linear_operator(cpr)
    solver = Jutul.Krylov.GmresSolver(op, rhs, 50)
    solve_fixed() = Jutul.Krylov.gmres!(solver, op, rhs;
        N = precond, rtol = 0.005, atol = 0.0, itmax = 100, history = true)
    linear_solve = term_row("fixed_system_gmres", solve_fixed; inner = 1)
    linear_solve["iterations"] = solver.stats.niter
    linear_solve["solved"] = solver.stats.solved
    mul!(x, op, solver.x, 1.0, 0.0)
    linear_solve["true_relative_residual"] = norm(rhs - x)/norm(rhs)
    push!(terms, linear_solve)
    H = cpr.pressure_precond.factor
    push!(terms, term_row("amg_pattern_check", () -> KP.same_pattern(H, H.levels[1].A)))
    for (i, level) in enumerate(H.levels)
        row = Dict("index" => i, "rows" => size(level.A, 1), "terms" => Any[])
        push!(row["terms"], term_row("update$(i)_smoother", () ->
            KP.update_level_smoother!(level.smoother, level.A, H.options)))
        if i < length(H.levels)
            coarse, P, G = H.levels[i + 1].A, level.P, level.galerkin
            row["galerkin_triples"] = length(G.a_index)
            row["galerkin_bytes"] = sum(sizeof(getproperty(G, f)) for f in fieldnames(typeof(G)))
            push!(row["terms"], term_row("update$(i)_galerkin", () ->
                KP.galerkin!(coarse, level.A, P, G)))
        elseif !isnothing(level.coarse_solver)
            push!(row["terms"], term_row("coarse_solver_update", () ->
                KP.update_coarse_solver!(level.coarse_solver, level.A)))
        end
        push!(updates, row)
    end
    result = Dict("fixture" => ARGS[1], "threads" => Threads.nthreads(:default),
        "julia" => string(VERSION), "pinned_cpus" => get(ENV, "PERF_PIN_CPUS", ""),
        "factor_comparison" => "native and level ILU are both refactored from the saved assembly matrix",
        "terms" => terms, "schur" => couplings, "amg_update" => updates,
        "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`),
        "darcy_commit" => readchomp(`git -C $(pkgdir(JutulDarcy)) rev-parse HEAD`))
    open(ARGS[2], "w") do io
        TOML.print(io, result)
    end
end
probe_full_terms()
