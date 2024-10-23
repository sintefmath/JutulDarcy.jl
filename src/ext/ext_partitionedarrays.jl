"""
    simulate_reservoir_parray(case, mode = :mpi; kwarg...)

Run simulation with parray. This function is primarily for testing.
[`simulate_reservoir`](@ref) can do the same job by passing the correct mode.
"""
function simulate_reservoir_parray(case, mode = :mpi; kwarg...)
    sim, cfg = setup_reservoir_simulator(case; mode = mode, kwarg...)
    return simulate!(sim, case.dt, forces = case.forces, config = cfg)
end

function setup_reservoir_simulator_parray

end
