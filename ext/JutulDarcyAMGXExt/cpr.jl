function JutulDarcy.update_amgx_pressure_system!(amgx::AMGXPreconditioner, A_p::Jutul.StaticCSR.StaticSparsityMatrixCSR)
    data = amgx.data
    initialize_amgx_structure!(amgx, A_p)
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

function initialize_amgx_structure!(amgx::AMGXPreconditioner, A::Jutul.StaticCSR.StaticSparsityMatrixCSR{Tv, <:Any}) where Tv
    if !haskey(amgx.data, :storage)
        if Tv == Float64
            amgx_mode = AMGX.dDDI
        else
            amgx_mode = AMGX.dFFI
        end
        # TODO: Tv isn't really coming from the matrix, it should be set beforehand in the krylov solver.
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
    end
    return amgx
end

