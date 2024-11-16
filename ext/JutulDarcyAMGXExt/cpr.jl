function JutulDarcy.update_amgx_pressure_system!(amgx::AMGXPreconditioner, A::Jutul.StaticCSR.StaticSparsityMatrixCSR, Tv, cpr::JutulDarcy.CPRPreconditioner, recorder)
    do_p_update = JutulDarcy.should_update_cpr(cpr, recorder, :amg)
    already_setup = haskey(amgx.data, :storage)
    if already_setup
        if do_p_update || cpr.partial_update
            pb = amgx.data[:buffer_p]
            s = amgx.data[:storage]
            A_gpu = s.matrix
            @assert nnz(A) == nnz(A_gpu)
            @tic "AMGX coefficients" AMGX.replace_coefficients!(A_gpu, A.At.nzval)
        end
    else
        if Tv == Float64
            amgx_mode = AMGX.dDDI
        else
            amgx_mode = AMGX.dFFI
        end
        n, m = size(A)
        @assert n == m
        s = AMGXStorage(amgx.settings, amgx_mode)
        # RHS and solution vectors to right size just in case
        AMGX.set_zero!(s.x, n)
        AMGX.set_zero!(s.r, n)

        row_ptr = Cint.(A.At.colptr .- 1)
        colval = Cint.(A.At.rowval .- 1)
        nzval = A.At.nzval
        # TODO: Support for other types than Float64, should be done in setup of
        # pressure system
        AMGX.pin_memory(nzval)
        AMGX.upload!(s.matrix,
            row_ptr,
            colval,
            nzval
        )
        amgx.data[:nzval] = A.At.nzval
        amgx.data[:storage] = s
        amgx.data[:block_size] = 1
        amgx.data[:n] = n
        amgx.data[:buffer_p] = AMGX.CUDA.CuVector{Tv}(undef, n)
    end
    if do_p_update
        @tic "AMGX setup" if already_setup && amgx.resetup
            AMGX.resetup!(s.solver, s.matrix)
        else
            AMGX.setup!(s.solver, s.matrix)
        end
    end
    return amgx
end

mutable struct AMGXStorage{C, R, V, M, S}
    config::C
    resources::R
    r::V
    x::V
    matrix::M
    solver::S
    function AMGXStorage(settings::AbstractDict, amgx_mode = AMGX.dDDI)
        config = AMGX.Config(settings)
        resources = AMGX.Resources(config)
        r = AMGX.AMGXVector(resources, amgx_mode)
        x = AMGX.AMGXVector(resources, amgx_mode)
        matrix = AMGX.AMGXMatrix(resources, amgx_mode)
        solver = AMGX.Solver(resources, amgx_mode, config)
        function finalize_storage!(amgx_s::AMGXStorage)
            close(amgx_s.solver)
            close(amgx_s.x)
            close(amgx_s.r)
            close(amgx_s.matrix)
            close(amgx_s.resources)
            close(amgx_s.config)
        end
        s = new{
            typeof(config),
            typeof(resources),
            typeof(r),
            typeof(matrix),
            typeof(solver)
        }(config, resources, r, x, matrix, solver)
        return finalizer(finalize_storage!, s)
    end
end

function JutulDarcy.gpu_amgx_solve!(amgx::AMGXPreconditioner, r_p)
    Jutul.apply!(r_p, amgx, r_p)
end

function JutulDarcy.gpu_cpr_setup_buffers!(cpr, J_bsr, r_cu, op, recorder)
    # TODO: Rename this function
    data = cpr.pressure_precond.data
    is_first = !haskey(data, :w_p)
    cpr_s = cpr.storage
    if is_first
        data[:w_p_cpu] = CUDA.pin(cpr_s.w_p)
        Tv = eltype(r_cu)
        n = length(r_cu)
        bz = cpr_s.block_size
        cpr.pressure_precond.data[:buffer_full] = similar(r_cu)
        bz_w, n_w = size(cpr_s.w_p)
        # @assert n == bz*n_w
        data[:w_p] = AMGX.CUDA.CuMatrix{Tv}(undef, bz_w, n_w)
        data[:main_system] = J_bsr
        data[:operator] = op
    end
    do_p_update = JutulDarcy.should_update_cpr(cpr, recorder, :amg)
    if is_first || do_p_update || cpr.partial_update
        # Put updated weights on GPU
        w_p_cpu = data[:w_p_cpu]
        w_p_gpu = data[:w_p]
        @assert size(w_p_cpu) == size(w_p_gpu)
        @tic "weights to gpu" copyto!(w_p_gpu, w_p_cpu)
    end
    return cpr
end
