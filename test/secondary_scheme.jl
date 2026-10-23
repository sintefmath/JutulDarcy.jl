using Jutul, JutulDarcy, Test

@testset "Reservoir secondary-variable scheduling" begin
    domain = reservoir_domain(CartesianMesh((4, 1, 1), (4.0, 1.0, 1.0)))
    system = ImmiscibleSystem((LiquidPhase(), VaporPhase());
        reference_densities = [1000.0, 100.0])
    wells = [setup_vertical_well(domain, 1, 1; name = :Producer, simple_well = true)]
    for backend in (:csr, :csc), block_backend in (false, true)
        @testset "$backend, block_backend=$block_backend" begin
            model = setup_reservoir_model(domain, system; wells, backend, block_backend)
            @test secondary_scheme(model) == :flat
            @test all(m -> secondary_scheme(m) == :flat, values(model.models))
            for scheme in (:flat, :fused, :strided)
                model = setup_reservoir_model(domain, system;
                    wells, backend, block_backend, secondary_scheme = scheme,
                    minbatch = 2, nthreads = min(2, Threads.nthreads()))
                @test secondary_scheme(model) == scheme
                @test all(m -> secondary_scheme(m) == scheme, values(model.models))
                @test Jutul.is_cell_major(Jutul.matrix_layout(reservoir_model(model).context)) == block_backend
                @test Jutul.minbatch(reservoir_model(model).context) == 2
            end
            model = setup_reservoir_model(domain, system;
                wells, backend, block_backend, context = nothing, secondary_scheme = :fused)
            @test all(m -> secondary_scheme(m) == :fused, values(model.models))
        end
    end

    context = DefaultContext(secondary_scheme = :strided)
    reservoir_context = ParallelCSRContext(secondary_scheme = :fused)
    model = setup_reservoir_model(domain, system; wells, context, reservoir_context)
    @test reservoir_model(model).context === reservoir_context
    @test model[:Producer].context === context
    @test model[:Facility].context === context
    @test_throws ArgumentError setup_reservoir_model(domain, system; secondary_scheme = :unknown)
end

@testset "Secondary-variable scheduling in KA simulations" begin
    reference = nothing
    for scheme in (:flat, :fused, :strided)
        case = JutulDarcy.setup_mini_wellcase(:immiscible_2ph;
            nstep = 1, total_time = 0.1, secondary_scheme = scheme, minbatch = 1)
        simulator, config = setup_reservoir_simulator(case;
            mode = :ka, ka_backend = Jutul.KernelExecution.KernelAbstractions.CPU(),
            linear_solver = nothing, timesteps = :none, info_level = -1,
            error_on_incomplete = true)
        @test secondary_scheme(simulator.model) == scheme
        @test all(m -> secondary_scheme(m) == scheme, values(simulator.model.models))
        result = simulate_reservoir(case; simulator, config)
        @test length(result.states) == 1
        state = only(result.states)
        if isnothing(reference)
            reference = state
        else
            @test state[:Pressure] ≈ reference[:Pressure]
            @test state[:Saturations] ≈ reference[:Saturations]
        end
    end
end
