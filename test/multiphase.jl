using JutulDarcy, Jutul
using Test

function test_multiphase(grid = CartesianMesh((2, 2), (2.0, 2.0)); setup = "two_phase_simple", debug_level = 1, linear_solver = nothing, kwarg...)
    state0, model, prm, f, t = get_test_setup(grid, case_name = setup; kwarg...)
    sim = Simulator(model, state0 = state0, parameters = prm)
    if linear_solver != :auto
        arg = (linear_solver = linear_solver, )
    else
        arg = NamedTuple()
    end
    cfg = simulator_config(sim; info_level = -1, debug_level = debug_level, arg...)
    simulate(sim, t, forces = f, config = cfg)
    return true
end

bctx = DefaultContext(matrix_layout = BlockMajorLayout())
ctx = DefaultContext(matrix_layout = EntityMajorLayout())

setups = ["two_phase_simple",
          "two_phase_fake_wells",
          "three_phase_fake_wells",
          "spe1_like_fake_wells",
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
            @testset "Block assembly, ILU0" begin
                @test test_multiphase(setup = setup, context = bctx, linear_solver = GenericKrylov(preconditioner = ILUZeroPreconditioner()))
            end
            @testset "Block assembly, CPR" begin
                for strategy in [:quasi_impes, :true_impes]
                    prec = CPRPreconditioner(strategy = strategy)
                    @test test_multiphase(setup = setup, context = bctx, linear_solver = GenericKrylov(preconditioner = prec))
                end
            end
            @testset "Block assembly, auto" begin
                @test test_multiphase(setup = setup, context = bctx, linear_solver = :auto)
            end
            @testset "Unit major assembly, direct solver" begin
                @test test_multiphase(setup = setup, context = ctx)
            end
        end
    end
end
