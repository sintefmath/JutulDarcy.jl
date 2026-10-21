# Included into KAPreconditioners for an isolated CPU experiment. All smoothing
# and update kernels remain the existing SPAI0/GaussSeidel implementations.
struct PerfHybridSmoother <: AbstractSmoother
    steps::Int
    damping::Float64
    cutoff::Int
end
PerfHybridSmoother(cutoff = 16_384; steps = 1, damping = 1.0) =
    PerfHybridSmoother(steps, damping, cutoff)

perf_smoother_config(state::SPAI0State, config) = SPAI0State(
    state.diagonal, state.temporary, state.residual, config,
    state.backend, state.block_size, state.n)
perf_smoother_config(state::GaussSeidelState, config) = GaussSeidelState(
    state.matrix, state.inverse_diagonal, state.correction, state.residual,
    config, state.backend, state.block_size, state.n)

function setup_smoother(A::StaticSparsityMatrixCSR, config::PerfHybridSmoother;
        reuse = nothing, reallocation_tracker = nothing)
    matrix_backend(A) isa KernelAbstractions.CPU || error("CPU-only experiment")
    fine = matrix_nrows(A) > config.cutoff
    base = fine ? SPAI0(config.steps, config.damping) : GaussSeidel(config.steps, config.damping)
    compatible = fine ? reuse isa SPAI0State : reuse isa GaussSeidelState
    # Restore the concrete base config before invoking its existing reuse path,
    # then retain the hybrid config so ordinary numeric resets remain valid.
    old = compatible ? perf_smoother_config(reuse, base) : nothing
    state = setup_smoother(A, base; reuse = old, reallocation_tracker)
    return perf_smoother_config(state, config)
end
