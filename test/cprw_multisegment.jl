using Jutul, JutulDarcy, Test, JLArrays, CUDA
using Jutul.KernelExecution: KernelAbstractions

function cprw_multisegment_case(; merged = false, split_wells = false, thermal = false)
    bar = si_unit(:bar)
    domain = reservoir_domain(CartesianMesh((3, 1, 3), (300.0, 100.0, 45.0));
        permeability = 1e-13, porosity = 0.2)
    wells = Any[setup_vertical_well(domain, i, 1; name, simple_well = false)
        for (i, name) in enumerate((:Producer, :Injector))]
    push!(wells, setup_well(domain, [3]; name = :Simple, simple_well = true))
    system = ImmiscibleSystem((LiquidPhase(), VaporPhase());
        reference_densities = [1000.0, 100.0])
    model = setup_reservoir_model(domain, system; wells, split_wells, thermal, backend = :csr)
    state_options = thermal ? (Temperature = 350.0,) : NamedTuple()
    state0 = setup_reservoir_state(model; Pressure = 200bar, Saturations = [0.8, 0.2], state_options...)
    controls = Dict(
        :Producer => ProducerControl(BottomHolePressureTarget(190bar)),
        :Injector => InjectorControl(TotalRateTarget(1e-4), [1.0, 0.0]; density = 1000.0, temperature = 350.0),
        :Simple => DisabledControl(),
    )
    next_controls = copy(controls)
    next_controls[:Producer] = ProducerControl(TotalRateTarget(-1e-4))
    next_controls[:Injector] = InjectorControl(BottomHolePressureTarget(210bar), [1.0, 0.0]; density = 1000.0, temperature = 350.0)
    forces = [setup_reservoir_forces(model; control) for control in (controls, next_controls)]
    case = JutulCase(model, fill(0.05*si_unit(:day), 2), forces; state0)
    return merged ? JutulDarcy.merge_similar_wells(case) : case
end

@testset "CPRW multisegment wells" begin
    options = (info_level = -1, timesteps = :none, error_on_incomplete = true,
        failure_cuts_timestep = false, tol_cnv = 1e-7, tol_mb = 1e-9, tol_cnv_well = 1e-7)
    policies = [
        (mode = :default,),
        (mode = :ka, ka_backend = KernelAbstractions.CPU()),
        (mode = :ka, ka_backend = JLBackend()),
    ]
    if CUDA.functional()
        CUDA.allowscalar(false)
        push!(policies, (mode = :ka_cuda,))
    end
    # Exercise backend reductions on ordinary and merged wells. Split groups,
    # thermal unknowns and the other fluid systems each need one focused case.
    cases = [
        ("ordinary", cprw_multisegment_case(), true, 3),
        ("merged", cprw_multisegment_case(; merged = true), true, 3),
        ("split", cprw_multisegment_case(; split_wells = true), false, 3),
        ("thermal", cprw_multisegment_case(; split_wells = true, thermal = true), false, 3),
    ]
    for physics in (:compositional_2ph_3c, :bo_spe1)
        case = JutulDarcy.setup_mini_wellcase(Val(physics);
            simple_well = false, nstep = 1, total_time = si_unit(:day), backend = :csr)
        push!(cases, (string(physics), case, false, 2))
    end
    for (name, case, test_backends, nwell) in cases
        reference = simulate_reservoir(case; precond = :ilu0, options...)
        for policy in (test_backends ? policies : policies[1:1])
            @testset "$name, $policy" begin
                result = simulate_reservoir(case; precond = :cprw, policy..., options...)
                @test length(result.states) == length(reference.states) == length(case.dt)
                for (state, ref) in zip(result.states, reference.states)
                    for field in (:Pressure, :Saturations, :Temperature, :OverallMoleFractions, :ImmiscibleSaturation)
                        if haskey(ref, field)
                            @test state[field] ≈ ref[field] rtol = 1e-7 atol = 1e-10
                        end
                    end
                end
                for well in well_symbols(case.model), quantity in (:bhp, :rate)
                    @test result.wells[well, quantity] ≈ reference.wells[well, quantity] rtol = 1e-7 atol = 1e-10
                end
                # Ensure CPRW actually enriches and couples the pressure solve.
                cpr = result.extra[:config][:linear_solver].preconditioner
                pressure = Matrix(JutulDarcy.cpr_host_matrix(cpr.storage.A_p))
                ncell = number_of_cells(reservoir_model(case.model).domain)
                @test size(pressure) == (ncell + nwell, ncell + nwell)
                @test any(!iszero, pressure[1:ncell, (ncell + 1):end])
                @test any(!iszero, pressure[(ncell + 1):end, 1:ncell])
            end
        end
    end
end
