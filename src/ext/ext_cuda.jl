mutable struct CUDAReservoirKrylov{Tv, Ti} <: Jutul.AbstractKrylov
    solver::Symbol
    config::IterativeSolverConfig
    preconditioner
    data
    storage
end

function CUDAReservoirKrylov(solver = :gmres, prec = ILUZeroPreconditioner();
        Float_t = Float64,
        Int_t = Int32,
        kwarg...
    )
    cfg = IterativeSolverConfig(; kwarg...)
    CUDAReservoirKrylov{Float_t, Int_t}(solver, cfg, prec, Dict{Symbol, Any}(), nothing)
end

function Jutul.linear_solve!(lsys::Jutul.LSystem,
        krylov::CUDAReservoirKrylov{Tv, Ti},
        model,
        storage = nothing,
        dt = nothing,
        recorder = ProgressRecorder(),
        executor = default_executor();
        dx = lsys.dx_buffer,
        r = Jutul.vector_residual(lsys),
        atol = Jutul.linear_solver_tolerance(krylov, :absolute),
        rtol = Jutul.linear_solver_tolerance(krylov, :relative)
    ) where {Tv, Ti}
    cfg = krylov.config

    atol = convert(Tv, atol)
    rtol = convert(Tv, rtol)
    t_prep0 = @elapsed Jutul.prepare_linear_solve!(lsys)

    is_single_linear_system = lsys isa LinearizedSystem

    L = lsys[1,1]
    J = L.jac
    J::Jutul.StaticSparsityMatrixCSR
    csr_block_buffer = L.jac_buffer

    bz = Jutul.block_size(L)
    sz = size(J)
    n, m = sz
    @assert n == m
    is_first = !haskey(krylov.data, :J)
    t_setup = @elapsed if is_first
        krylov.data[:J], krylov.data[:r] = build_gpu_block_system(Ti, Tv, sz, bz, J.At.colptr, J.At.rowval, csr_block_buffer, r)
        krylov.data[:schur] = build_gpu_schur_system(Ti, Tv, bz, lsys)
        krylov.data[:dx_cpu] = zeros(n*bz)
    end
    J_bsr = krylov.data[:J]
    r_cu = krylov.data[:r]
    schur = krylov.data[:schur]
    dx_cpu = krylov.data[:dx_cpu]

    t_gpu_update = @elapsed begin
        if !is_first
            update_gpu_block_residual!(r_cu, bz, r)
            update_gpu_block_system!(J_bsr, bz, csr_block_buffer)
            update_gpu_schur_system!(schur, lsys)
        end
    end
    t_prep = t_gpu_update + t_setup + t_prep0
    max_it = krylov.config.max_iterations

    Jutul.update_preconditioner!(krylov.preconditioner, J_bsr, r_cu, model.context, executor)
    prec_op = Jutul.preconditioner(krylov, lsys, model, storage, recorder, Tv)
    if cfg.precond_side == :right
        preconditioner_arg = (N = prec_op, )
    else
        preconditioner_arg = (M = prec_op, )
    end
    operator = gpu_system_linear_operator(J_bsr, schur, Tv)
    solve_f, F = Jutul.krylov_jl_solve_function(krylov, operator, r_cu)
    solve_f(F, operator, r_cu;
        itmax = max_it,
        preconditioner_arg...,
        verbose = 0,
        ldiv = false,
        history = true,
        rtol = rtol,
        atol = atol
    )
    dx_gpu, stats = (krylov.storage.x, krylov.storage.stats)

    copyto!(dx_cpu, dx_gpu)
    Jutul.update_dx_from_vector!(lsys, dx_cpu, dx = dx)
    # @. dx = -dx

    res = stats.residuals
    solved = stats.solved
    n = length(res) - 1

    t_prec = 0.0
    if prec_op isa Jutul.PrecondWrapper
        t_prec_op = prec_op.time
        t_prec_count = prec_op.count
    else
        t_prec_op = 0.0
        t_prec_count = 0
    end
    return Jutul.linear_solve_return(solved, n, stats, precond = t_prec_op, precond_count = t_prec_count, prepare = t_prec + t_prep)
end

function Base.show(io::IO, krylov::CUDAReservoirKrylov)
    rtol = linear_solver_tolerance(krylov.config, :relative)
    atol = linear_solver_tolerance(krylov.config, :absolute)
    print(io, "CUDAReservoirKrylov using $(krylov.solver) (ϵₐ=$atol, ϵ=$rtol) with prec = $(typeof(krylov.preconditioner))")
end

function gpu_system_linear_operator(J, schur::Nothing, Tv)
    return Jutul.LinearOperators.LinearOperator(J)
end

function gpu_system_linear_operator(J, schur, Tv)
    n, m = size(J)
    C = schur[:C]
    D = schur[:D]
    E_factor = schur[:E_factor]
    buf = schur[:buf]
    buf_1 = schur[:buf_1_cpu]
    buf_2 = schur[:buf_2_cpu]

    mul!(x, y, α, β) = schur_mul_gpu!(x, y, α, β, J, C, D, buf, buf_1, buf_2, E_factor)
    return Jutul.LinearOperators.LinearOperator(Tv, n, m, false, false, mul!)
end

function schur_mul_gpu!

end

function build_gpu_block_system

end

function update_gpu_block_system!

end

function build_gpu_schur_system

end

function update_gpu_block_residual!

end

function update_gpu_schur_system!

end
