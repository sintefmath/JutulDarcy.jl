function JutulDarcy.update_amgx_pressure_system!(amgx::AMGXPreconditioner, A_p::Jutul.StaticCSR.StaticSparsityMatrixCSR, Tv)
    data = amgx.data
    initialize_amgx_structure!(amgx, A_p, Tv)
    update_pressure_system!(amgx, A_p)
    return amgx
end

mutable struct AMGXStorage{C, R, V, M, S}
    config::C
    resources::R
    r::V
    x::V
    matrix::M
    solver::S
    function AMGXStorage(config, amgx_mode = AMGX.dDDI)
        resources = AMGX.Resources(config)
        r = AMGX.AMGXVector(resources, amgx_mode)
        x = AMGX.AMGXVector(resources, amgx_mode)
        matrix = AMGX.AMGXMatrix(resources, amgx_mode)
        solver = AMGX.Solver(resources, amgx_mode, config)
        function finalize_storage!(amgx_s::AMGXStorage)
            close(amgx_s.solver)
            close(amgx_s.matrix)
            close(amgx_s.r)
            close(amgx_s.x)
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

function update_pressure_system!(amgx, A_cpu)
    A_gpu = amgx.data[:storage].matrix
    AMGX.replace_coefficients!(A_gpu, A_cpu.At.nzval)
    return amgx
end

function initialize_amgx_structure!(amgx::AMGXPreconditioner, A::Jutul.StaticCSR.StaticSparsityMatrixCSR, Tv)
    if !haskey(amgx.data, :storage)
        if Tv == Float64
            amgx_mode = AMGX.dDDI
        else
            amgx_mode = AMGX.dFFI
        end
        n, m = size(A)
        @assert n == m
        config = AMGX.Config(amgx.settings)
        s = AMGXStorage(config, amgx_mode)
        # RHS and solution vectors to right size just in case
        AMGX.set_zero!(s.x, n)
        AMGX.set_zero!(s.r, n)

        row_ptr = Cint.(A.At.colptr .- 1)
        colval = Cint.(A.At.rowval .- 1)
        nzval = A.At.nzval
        nzval = Tv.(nzval)

        AMGX.upload!(s.matrix, 
            row_ptr,
            colval,
            nzval
        )
        amgx.data[:storage] = s
        amgx.data[:buffer_p] = AMGX.CUDA.CuVector{Tv}(undef, n)
    end
    return amgx
end

function JutulDarcy.gpu_cpr_setup_buffers!(cpr, J_bsr, r_cu, op)
    # TODO: Rename this.
    data = cpr.pressure_precond.data
    if !haskey(data, :w_p)
        Tv = eltype(r_cu)
        n = length(r_cu)
        cpr_s = cpr.storage
        bz = cpr_s.block_size

        # cpr.pressure_precond.data[:buffer_full] = AMGX.CUDA.CuVector{Tv}(undef, n)
        bz_w, n_w = size(cpr_s.w_p)
        @assert n == bz*n_w
        data[:w_p] = AMGX.CUDA.CuMatrix{Tv}(undef, bz_w, n_w)
        data[:main_system] = J_bsr
        data[:operator] = op
    end
    # Put updated weights on GPU
    w_p_cpu = cpr_s.w_p
    w_p_gpu = data[:w_p]
    @assert size(w_p_cpu) == size(w_p_gpu)
    copyto!(w_p_gpu, w_p_cpu)
    return cpr
end
