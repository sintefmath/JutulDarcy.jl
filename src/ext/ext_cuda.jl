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

function Jutul.linear_solve!(lsys::Jutul.LinearizedSystem,
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

    L = lsys[1,1]
    J = L.jac
    J::Jutul.StaticSparsityMatrixCSR
    r0 = L.r_buffer
    bz = Jutul.block_size(lsys)
    sz = size(J)
    n, m = sz
    @assert n == m
    is_first = !haskey(krylov.data, :J)
    t_setup = @elapsed if is_first
        krylov.data[:J], krylov.data[:r] = build_gpu_block_system(Ti, Tv, sz, bz, J.At.colptr, J.At.rowval, L.jac_buffer, r0)
    end
    J_bsr = krylov.data[:J]
    r_cu = krylov.data[:r]

    t_gpu_update = @elapsed begin
        update_gpu_block_residual!(r_cu, bz, r0)
        update_gpu_block_system!(J_bsr, bz, L.jac_buffer)
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
    solve_f, F = Jutul.krylov_jl_solve_function(krylov, J_bsr, r_cu)
    solve_f(F, J_bsr, r_cu;
        itmax = max_it,
        preconditioner_arg...,
        verbose = 0,
        ldiv = false,
        history = true,
        rtol = rtol,
        atol = atol
    )
    dx_gpu, stats = (krylov.storage.x, krylov.storage.stats)

    copyto!(dx, dx_gpu)
    @. dx = -dx

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

function build_gpu_block_system

end

function update_gpu_block_system!

end

function update_gpu_block_residual!

end
