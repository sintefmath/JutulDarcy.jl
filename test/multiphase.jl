using JutulDarcy, Jutul
using Test

function test_multiphase(grid = "pico"; setup = "two_phase_simple", debug_level = 1, linear_solver = nothing, kwarg...)
    state0, model, prm, f, t = get_test_setup(grid, case_name = setup; kwarg...)
    sim = Simulator(model, state0 = state0, parameters = prm)
    if linear_solver == :auto
        linear_solver = reservoir_linsolve(model)
    end
    cfg = simulator_config(sim, info_level = -1, debug_level = debug_level, linear_solver = linear_solver)
    simulate(sim, t, forces = f, config = cfg)
    return true
end

bctx = DefaultContext(matrix_layout = BlockMajorLayout())
ctx = DefaultContext(matrix_layout = UnitMajorLayout())

setups = ["two_phase_simple",
          "two_phase_fake_wells",
          "three_phase_fake_wells",
          "simple_compositional_fake_wells",
          "compositional_three_phases"
          ]

@testset "Multi phase flow" begin
    for setup in setups
        @testset "$setup" begin
            @testset "Default, direct solver" begin
                @test test_multiphase(setup = setup)
            end
            @testset "Default backend, Krylov solver" begin
                @test test_multiphase(setup = setup, linear_solver = GenericKrylov())
            end
            @testset "Block assembly, Krylov solver" begin
                @test test_multiphase(setup = setup, context = bctx, linear_solver = GenericKrylov())
            end
            @testset "Block assembly, preconditioners" begin
                @test test_multiphase(setup = setup, context = bctx, linear_solver = GenericKrylov(preconditioner = ILUZeroPreconditioner()))
                @test test_multiphase(setup = setup, context = bctx, linear_solver = GenericKrylov(preconditioner = CPRPreconditioner()))
                @test test_multiphase(setup = setup, context = bctx, linear_solver = :auto)
            end
            @testset "Unit major assembly, direct solver" begin
                @test test_multiphase(setup = setup, context = ctx)
            end
        end
    end
end
