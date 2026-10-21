# Compare natural and saved dependency traversals on the same Float64 factors.
include("probe_terms.jl")

function serial_factor!(state, A, schedule)
    copyto!(state.factors, 1, A.nzval, 1, KP.matrix_nonzeros(A))
    return KP.ilu_factor_cpu!(state, schedule.lower_offsets, schedule.lower_rows)
end

function probe_serial_ilu()
    fixture = deserialize(ARGS[1])
    state = fixture.preconditioner.system_precond.factor
    A = JutulDarcy.reservoir_jacobian(fixture.linear_system)
    @assert eltype(eltype(A)) === Float64
    level_order = (lower_offsets = state.factor_offsets, lower_rows = state.factor_rows,
        upper_offsets = state.upper_offsets, upper_rows = state.upper_rows)
    selected = KP.cpu_smoother_schedule(state)
    serial_factor!(state, A, level_order)
    factors, diagonal = copy(state.factors), copy(state.inverse_diagonal)
    serial_factor!(state, A, selected)
    factor_error = norm(state.factors - factors)/norm(factors)
    diagonal_error = norm(state.inverse_diagonal - diagonal)/norm(diagonal)
    @assert factor_error < 1e-14 && diagonal_error < 1e-14
    rhs = copy(state.work)
    fill!(rhs, eltype(rhs)(ntuple(_ -> 1.0, length(first(rhs)))))
    x = similar(rhs)
    KP.ilu_solve_cpu!(x, state, rhs, state.work, level_order)
    reference = copy(x)
    KP.ilu_solve_cpu!(x, state, rhs, state.work)
    solve_error = norm(x - reference)/norm(reference)
    @assert solve_error < 1e-14
    terms = Any[]
    for repeat in 1:3
        choices = (("levels", level_order), ("selected", selected))
        for (label, schedule) in (isodd(repeat) ? choices : reverse(choices))
            push!(terms, merge(Dict("repeat" => repeat),
                term_row("$(label)_apply", () -> KP.ilu_solve_cpu!(x, state, rhs, state.work, schedule))))
            push!(terms, merge(Dict("repeat" => repeat),
                term_row("$(label)_factor", () -> serial_factor!(state, A, schedule); inner=3)))
        end
    end
    open(ARGS[2], "w") do io
        TOML.print(io, Dict("fixture" => ARGS[1], "threads" => Threads.nthreads(:default),
            "factor_relative_error" => factor_error, "diagonal_relative_error" => diagonal_error,
            "solve_relative_error" => solve_error, "terms" => terms,
            "jutul_commit" => readchomp(`git -C $(pkgdir(Jutul)) rev-parse HEAD`)))
    end
end
probe_serial_ilu()
