using Jutul
using JutulDarcy
using Test
using Adapt

const f32_backend = Jutul.KernelExecution.KernelAbstractions.CPU()
const f32_step = 0.01*si_unit(:day)

function short_reservoir_case(case)
    return JutulCase(case.model, [f32_step], case.forces[1:1];
        state0 = case.state0, parameters = case.parameters)
end

function test_immiscible_secondary_float32(simulator)
    model = simulator.model.models[:Reservoir]
    variables = model.secondary_variables
    original_state = Jutul.evaluation_state(simulator.storage.Reservoir)
    state = merge(Jutul.data(original_state),
        (StaticFluidVolume = original_state.FluidVolume,))
    pore_volume = Adapt.adapt(model.context,
        JutulDarcy.LinearlyCompressiblePoreVolume(
            reference_pressure = 1.0e5, expansion = 1.0e-9))

    @test model.system isa JutulDarcy.ImmiscibleSystem
    @test typeof(Jutul.value(first(state.Pressure))) === Float32
    @test typeof(first(JutulDarcy.reference_densities(model.system))) === Float32
    @test typeof(first(pore_volume.reference_pressure)) === Float32
    @test typeof(first(pore_volume.expansion)) === Float32
    @test typeof(first(variables.RelativePermeabilities.krw).connate) === Float32
    for phase_table in variables.PhaseViscosities.pvt
        @test typeof(first(phase_table.tab).p_ref) === Float32
    end

    nph = length(model.system.phases)
    for name in (:FluidVolume, :PhaseMassDensities, :TotalMasses,
            :PhaseViscosities, :RelativePermeabilities,
            :PhaseMobilities, :PhaseMassMobilities)
        variable = name === :FluidVolume ? pore_volume : variables[name]
        output = if name === :FluidVolume
            Vector{Any}(undef, 1)
        else
            Matrix{Any}(undef, nph, 1)
        end
        Jutul.update_secondary_variable!(output, variable, model, state, 1:1)
        values = if name === :FluidVolume
            (output[1],)
        else
            (output[ph, 1] for ph in 1:nph)
        end
        @test all(typeof(value) === eltype(state.Pressure) for value in values)
    end
end

function test_float32_reservoir_case(case, linear_float, linear_index;
        precond = :cpr, check_secondary_types = false)
    simulator, config = setup_reservoir_simulator(case;
        mode = :ka,
        ka_backend = f32_backend,
        float_type = Float32,
        index_type = Int32,
        linear_float_type = linear_float,
        linear_index_type = linear_index,
        precond = precond,
        info_level = -1,
        failure_cuts_timestep = false,
        error_on_incomplete = true)

    context = simulator.model.context
    system = simulator.storage.LinearizedSystem
    pressure = simulator.storage.Reservoir.state.Pressure
    @test Jutul.float_type(context) === Float32
    @test Jutul.index_type(context) === Int32
    @test Jutul.linear_float_type(context) === linear_float
    @test Jutul.linear_index_type(context) === linear_index
    @test typeof(Jutul.value(first(pressure))) === Float32
    @test eltype(system.r_buffer) === linear_float
    @test eltype(system[1, 1].jac.rowptr) === linear_index
    @test config[:tolerances][:Facility][:control_equation].Abs == 0.01

    if check_secondary_types
        test_immiscible_secondary_float32(simulator)
    end

    result = simulate_reservoir(case; simulator, config)
    solver = config[:linear_solver]
    @test eltype(solver.storage.x) === linear_float
    if solver.preconditioner isa CPRPreconditioner
        @test eltype(solver.preconditioner.storage.w_rhs) === linear_float
    end
    @test length(result.states) == 1
    if length(result.states) == 1
        @test all(isfinite, result.states[1][:Pressure])
    end
end

@testset "Float32 two-phase hysteresis" begin
    context = KernelAbstractionsContext(f32_backend;
        float_type = Float32, index_type = Int32)
    drain = PhaseRelativePermeability(
        Float32[0, 0.1, 0.3, 1], Float32[0, 0, 0.2, 1])
    imb = PhaseRelativePermeability(
        Float32[0, 0.2, 0.4, 1], Float32[0, 0, 0.2, 1])
    relperm = JutulDarcy.ReservoirRelativePermeabilities(
        w = (drain, imb), ow = (drain, imb),
        hysteresis_w = JutulDarcy.KilloughHysteresis())
    relperm = Adapt.adapt(context, relperm)
    system = ImmiscibleSystem(:wo;
        reference_densities = (1000.0f0, 800.0f0))
    saturation = Jutul.ForwardDiff.Dual{Nothing}(0.5f0, 1.0f0)
    state = (
        Saturations = reshape([saturation, one(saturation) - saturation], 2, 1),
        MaxSaturations = reshape(Float32[0.8, 0.8], 2, 1))
    output = Matrix{Any}(undef, 2, 1)
    Jutul.update_secondary_variable!(output, relperm, (; system), state, 1:1)
    @test all(typeof(output[phase, 1]) === typeof(saturation) for phase in 1:2)
    @test isfinite(Jutul.ForwardDiff.partials(output[1, 1])[1])
end

@testset "Float32 immiscible WATDENT density" begin
    context = KernelAbstractionsContext(f32_backend;
        float_type = Float32, index_type = Int32)
    water = JutulDarcy.PVTW(ConstMuBTable(
        Float32[1.0e5, 1, 1.0e-9, 1.0e-3, 0]))
    oil = JutulDarcy.PVCDO(
        [Float32[1.0e5, 1, 1.0e-9, 2.0e-3, 0]])
    @test typeof(first(water.tab).b_ref) === Float32
    @test typeof(first(oil.tab).b_ref) === Float32
    density = DeckPhaseMassDensities((water, oil);
        watdent = JutulDarcy.WATDENT([(300.0, 0.01, 0.001)]))
    density = Adapt.adapt(context, density)
    @test typeof(first(density.watdent.tab).T) === Float32
    system = ImmiscibleSystem(:wo;
        reference_densities = (1000.0f0, 800.0f0))
    pressure = Jutul.ForwardDiff.Dual{Nothing}(Float32(1.2e5), 1.0f0)
    temperature = Jutul.ForwardDiff.Dual{Nothing}(310.0f0, 1.0f0)
    output = Matrix{Any}(undef, 2, 1)
    JutulDarcy.update_deck_density!(output, density, (; system),
        [pressure], [temperature], 1:1)
    @test all(typeof(output[phase, 1]) === typeof(pressure) for phase in 1:2)
end

@testset "Float32 reservoir deck steps" begin
    for name in ("EGG", "SPE1", "SPE9")
        @testset "$name" begin
            path = JutulDarcy.GeoEnergyIO.test_input_file_path(
                name, "$name.DATA")
            case = setup_case_from_data_file(path;
                backend = :csr, block_backend = true)[1:1]
            short_case = short_reservoir_case(case)
            for (linear_float, linear_index) in
                    ((Float32, Int32), (Float64, Int64))
                @testset "linear $linear_float/$linear_index" begin
                    test_float32_reservoir_case(
                        short_case, linear_float, linear_index;
                        check_secondary_types =
                            name == "EGG" && linear_float === Float32)
                end
            end
        end
    end
end

@testset "Float32 mini reservoir steps" begin
    for physics in (:single_phase, :immiscible_2ph, :bo_spe1,
            :compositional_2ph_3c, :geothermal)
        @testset "$physics" begin
            case = JutulDarcy.setup_mini_wellcase(Val(physics);
                nstep = 1, total_time = f32_step,
                backend = :csr, block_backend = true)
            precond = :cpr
            if physics == :single_phase
                precond = :ka_spai0
            end
            for (linear_float, linear_index) in
                    ((Float32, Int32), (Float64, Int64))
                @testset "linear $linear_float/$linear_index" begin
                    test_float32_reservoir_case(
                        case, linear_float, linear_index;
                        precond = precond)
                end
            end
        end
    end
end
