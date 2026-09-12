"""
    setup_reservoir_linear_solver(model, precond = :cpr; <keyword arguments>)

Set up iterative linear solver for a reservoir model from [`setup_reservoir_model`](@ref).

# Arguments
- `model`: Reservoir model that will linearize the equations for the linear solver
- `precond=:cpr`: Preconditioner type to use. `:cpr` and `:cprw` select
  constrained-pressure-residual variants; smoother-only choices include
  `:ilu0`, `:jacobi`, `:spai0`, `:ka_ilu0`, `:ka_dilu`, and `:ka_spai0`.
- `backend=:auto`: Use the KernelAbstractions solver for a cell-major model
  with a KA context and the CPU/default fallback otherwise. `:ka`, `:cpu`, and
  the legacy transfer-based `:cuda` path can be selected explicitly.
- `v=0`: verbosity (can lead to a large amount of output)
- `solver=:bicgstab`: the symbol of a Krylov.jl solver (typically :gmres or :bicgstab)
- `update_interval=:once`: how often the CPR AMG hierarchy is reconstructed (:once, :iteration, :ministep, :step)
- `update_interval_partial=:iteration`: how often the pressure system is updated in CPR
- `max_coarse`: max size of coarse level if using AMG
- `amg_type`: pressure preconditioner implementation. `:ka` selects the new
  backend-portable AMG. An initialized Jutul preconditioner can also be passed.
- `smoother_type`: full-system smoother. The `:ka_*` choices use the new
  backend-portable smoothers.
- `amg_arg`, `smoother_arg`, `cpr_arg`: keyword arguments forwarded to the
  pressure AMG, full-system smoother, and CPR constructors, respectively.
- `cpr_type=nothing`: type of CPR (`:true_impes`, `:quasi_impes` or `nothing` for automatic)
- `partial_update=true`: perform partial update of CPR preconditioner outside of AMG update (see above)
- `rtol=1e-3`: relative tolerance for the linear solver
- `max_iterations=100`: limit for linear solver iterations

Additional keywords are passed onto the linear solver constructor.
"""
function setup_reservoir_linear_solver(model, precond = :cpr;
        backend = :auto,
        rtol = nothing,
        atol = nothing,
        v = 0,
        mode = :forward,
        solver = :bicgstab,
        max_iterations = nothing,
        update_interval = :iteration,
        update_interval_partial = :iteration,
        amg_type = nothing,
        smoother_type = :ilu0,
        max_coarse = 10,
        cpr_type = nothing,
        partial_update = update_interval == :once,
        amg_arg = NamedTuple(),
        smoother_arg = NamedTuple(),
        cpr_arg = NamedTuple(),
        precond_side = missing,
        float_type = Float64,
        kwarg...
    )
    is_equation_major = !Jutul.is_cell_major(matrix_layout(model.context))
    if backend == :auto
        backend = model.context isa Jutul.KernelAbstractionsContext &&
            !is_equation_major ?
            :ka : :cpu
    end
    backend in (:cpu, :cuda, :ka) || throw(ArgumentError(
        "Backend $backend not supported, must be :auto, :cpu, :ka or :cuda."))
    is_cpr = precond == :cpr || precond == :cprw
    if backend == :ka
        !is_equation_major || throw(ArgumentError(
            "Equation-major storage is not supported for KernelAbstractions solvers."))
        solver != :lu || throw(ArgumentError(
            "A direct LU solver is not supported for KernelAbstractions backends."))
        if isnothing(amg_type)
            amg_type = :ka
        end
        if smoother_type == :ilu0
            smoother_type = :ka_ilu0
        end
        krylov_constructor = GenericKrylov
        krylov_arg = NamedTuple()
    elseif backend == :cuda
        # Check assumptions
        !is_equation_major || throw(ArgumentError("Equation-major storage not supported for CUDA backend."))
        solver != :lu || throw(ArgumentError("LU direct solver not supported for CUDA backend."))
        has_cuda = !isnothing(Base.get_extension(JutulDarcy, :JutulDarcyCUDAExt))
        has_cuda || throw(ArgumentError("CUDA backend not available. You must run \"using CUDA\" before using this function."))
        has_amgx = !isnothing(Base.get_extension(JutulDarcy, :JutulDarcyAMGXExt))
        # Make sure that options are compatible with CUDA backend
        if is_cpr && !has_amgx
            jutul_message("AMGX", "AMGX not available, disabling CPR and falling back to ILU(0) preconditioner.")
            precond = :ilu0
        else
            amg_type = :amgx
        end
        if smoother_type != :ilu0
            jutul_message("CUDA", "Smoother $smoother_type not supported for CUDA, falling back to ILU(0).")
            smoother_type = :ilu0
        end
        krylov_constructor = CUDAReservoirKrylov
        krylov_arg = (Float_t = float_type, )
    else
        if solver == :lu
            return LUSolver()
        end
        if is_equation_major
            return nothing
        end
        krylov_constructor = GenericKrylov
        krylov_arg = NamedTuple()
        @assert float_type == Float64 "Only Float64 supported for CPU backend."
    end
    isnothing(amg_type) && (amg_type = default_amg_symbol())

    default_tol = 0.01
    max_it = 200
    if is_cpr
        if isnothing(cpr_type)
            if isa(model.system, ImmiscibleSystem)
                cpr_type = :analytical
            else
                cpr_type = :true_impes
            end
        end
        p_solve = default_psolve(; max_coarse = max_coarse, type = amg_type, amg_arg...)
        s = reservoir_system_preconditioner(smoother_type; smoother_arg...)
        prec = CPRPreconditioner(
            p_solve, s;
            strategy = cpr_type,
            variant = precond,
            update_interval = update_interval,
            partial_update = partial_update,
            update_interval_partial = update_interval_partial,
            mode = mode,
            cpr_arg...
        )
        default_tol = 0.005
        max_it = 50
    elseif precond in (:ilu0, :jacobi, :spai0, :ka_ilu0, :ka_dilu,
            :ka_spai0)
        selected = backend == :ka && precond == :ilu0 ? :ka_ilu0 : precond
        prec = reservoir_system_preconditioner(selected; smoother_arg...)
    else
        if precond isa Symbol
            error("Preconditioner $precond not supported for $(model.context)")
        else
            prec = precond
        end
    end
    if isnothing(rtol)
        rtol = default_tol
    end
    if isnothing(max_iterations)
        max_iterations = max_it
        if mode == :adjoint
            # No outer loop to control - add more iterations.
            max_iterations *= 4
        end
    end
    if ismissing(precond_side)
        if mode == :forward
            precond_side = :right
        else
            precond_side = :left
        end
    end
    lsolve = krylov_constructor(
        solver;
        verbose = v,
        preconditioner = prec,
        relative_tolerance = rtol,
        absolute_tolerance = atol,
        max_iterations = max_iterations,
        precond_side = precond_side,
        krylov_arg...,
        kwarg...
    )
    return lsolve
end

function reservoir_system_preconditioner(type; kwarg...)
    if type == :ilu0
        return ILUZeroPreconditioner(; kwarg...)
    elseif type == :jacobi
        return JacobiPreconditioner(; kwarg...)
    elseif type == :spai0
        return SPAI0Preconditioner(; kwarg...)
    elseif type == :ka_ilu0
        return Jutul.KASmootherPreconditioner(:ilu0; kwarg...)
    elseif type == :ka_dilu
        return Jutul.KASmootherPreconditioner(:dilu; kwarg...)
    elseif type == :ka_spai0
        return Jutul.KASmootherPreconditioner(:spai0; kwarg...)
    elseif type isa JutulPreconditioner
        return type
    else
        throw(ArgumentError("Unsupported reservoir smoother: $type"))
    end
end

function setup_reservoir_linear_solver(model::MultiModel, arg...; kwarg...)
    rmodel = reservoir_model(model)
    return setup_reservoir_linear_solver(rmodel, arg...; kwarg...)
end

function default_amg_symbol()
    if Jutul.check_hypre_availability(throw = false)
        amg_type = :hypre
    else
        amg_type = :smoothed_aggregation
    end
    return amg_type
end


function Jutul.select_linear_solver(m::SimulationModel{<:Any, S, <:Any, <:Any}; kwarg...) where S<:MultiPhaseSystem
    return setup_reservoir_linear_solver(m; kwarg...)
end
