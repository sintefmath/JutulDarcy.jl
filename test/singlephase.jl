using Jutul, JutulDarcy
using Test

function test_single_phase(grid = CartesianMesh((2, 2), (2.0, 2.0)); linear_solver = nothing, kwarg...)
    state0, model, prm, f, t = get_test_setup(grid, case_name = "single_phase_simple"; kwarg...)
    sim = Simulator(model, state0 = state0, parameters = prm)
    cfg = simulator_config(sim, info_level = -1)
    cfg[:linear_solver] = linear_solver
    simulate(sim, t, forces = f, config = cfg)
    # We just return true. The test at the moment just makes sure that the simulation runs.
    return true
end
@testset "Single-phase" begin
    @test test_single_phase()
end

agg = AMGPreconditioner(:smoothed_aggregation)
rs = AMGPreconditioner(:ruge_stuben)
@testset "Single-phase linear solvers" begin
    @test test_single_phase(linear_solver = GenericKrylov(preconditioner = agg))
    @test test_single_phase(linear_solver = GenericKrylov(preconditioner = rs))
    @test test_single_phase(linear_solver = GenericKrylov())
    @test test_single_phase(linear_solver = GenericKrylov(preconditioner = ILUZeroPreconditioner()))
end
