using Test, Jutul, JutulDarcy, HYPRE, PartitionedArrays
# Test distributed parallel solve using PartitionedArrays + HYPRE.
function setup_well_case(nx = 5, backend = :csr)
    bar = 1e5
    day = 3600*24.0
    Darcy = 9.869232667160130e-13

    ny = nx
    nz = 1
    dims = (nx, ny, nz)
    g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))

    domain = reservoir_domain(g, permeability = 0.1*Darcy, porosity = 0.2)
    Prod = setup_vertical_well(domain, 1, 1, name = :Producer);
    Inj = setup_well(domain, [(nx, ny, 1)], name = :Injector);
    phases = (LiquidPhase(), VaporPhase())
    rhoLS = 1000.0
    rhoGS = 100.0
    rhoS = [rhoLS, rhoGS]
    sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    model, parameters = setup_reservoir_model(
        domain, sys,
        wells = [Inj, Prod],
        backend = backend,
        split_wells = true
        )
    c = [1e-6/bar, 1e-4/bar]
    ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    replace_variables!(model, PhaseMassDensities = ρ);
    state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])
    dt = repeat([7.5]*day, 12*20)
    pv = pore_volume(model, parameters)
    inj_rate = sum(pv)/sum(dt)

    rate_target = TotalRateTarget(inj_rate)
    I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = rhoGS)
    bhp_target = BottomHolePressureTarget(50*bar)
    P_ctrl = ProducerControl(bhp_target)
    controls = Dict()
    controls[:Injector] = I_ctrl
    controls[:Producer] = P_ctrl

    forces = setup_reservoir_forces(model, control = controls)
    return JutulCase(model, dt, forces, state0 = state0, parameters = parameters)
end

function setup_bl_case(nc, backend = :csr; nstep = 10*nc)
    time = 1.0

    T = time
    tstep = repeat([T/nstep], nstep)
    domain = get_1d_reservoir(nc)
    nc = number_of_cells(domain)
    timesteps = tstep*3600*24
    bar = 1e5
    p0 = 100*bar
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()))
    if backend == :csr
        ctx = ParallelCSRContext(1, matrix_layout = BlockMajorLayout())
    else
        ctx = DefaultContext(matrix_layout = BlockMajorLayout())
    end
    model = SimulationModel(domain, sys, context = ctx)
    kr = BrooksCoreyRelativePermeabilities(sys, [2.0, 2.0], [0.2, 0.2])
    replace_variables!(model, RelativePermeabilities = kr)
    tot_time = sum(timesteps)
    pv = pore_volume(domain)
    irate = 500*sum(pv)/tot_time
    src  = SourceTerm(1, irate, fractional_flow = [1.0, 0.0])
    bc = FlowBoundaryCondition(nc, p0/2)
    forces = setup_forces(model, sources = src, bc = bc)
    parameters = setup_parameters(model, PhaseViscosities = [1e-3, 5e-3]) # 1 and 5 cP
    state0 = setup_state(model, Pressure = p0, Saturations = [0.0, 1.0])
    return JutulCase(model, tstep, forces, state0 = state0, parameters = parameters)
end
##
function compare_states(states_ref, states; rtol = 1e-4)
    for (state, state_ref) in zip(states, states_ref)
        for (k, v) in state_ref
            if v isa AbstractArray && eltype(v)<:Number
                @test v ≈ state[k] rtol = rtol
            end
        end
    end
end
##
function compare_ws(ws, ws_ref)
    @test keys(ws_ref) == keys(ws)
    for (k, v) in ws_ref
        if v isa AbstractArray && eltype(v)<:Number
            @test v ≈ ws[k]
        end
    end
end
##
@testset "PArray" begin
    num_procs_to_test = 1:5
    @testset "SimulationModel" begin
        for backend in [:csc, :csr]
            @testset "$backend" begin
                case = setup_bl_case(50, backend)
                arg = (info_level = -1,)
                states, reports = simulate(case; arg...)
                @testset "PArray native" begin
                    for np = num_procs_to_test
                        for order in [:default]#, :symrcm]
                            states_p, reports_p = simulate_reservoir_parray(case, 
                                :parray;
                                arg...,
                                parray_arg = (np = np, order = order),
                                output_path = tempdir(),
                                precond = :ilu0,
                                rtol = 1e-4,
                                timesteps = :none
                                )
                            compare_states(states, states_p, rtol = 1e-4)
                        end
                    end
                end
                @testset "PArray MPI" begin
                    states_m, reports_m = simulate_reservoir_parray(case; mode = :mpi, timesteps = :none, arg...)
                    compare_states(states, states_m)
                end
            end
        end
    end
    @testset "MultiModel" begin
        for backend in [:csc, :csr]
            @testset "$backend" begin
                case = setup_well_case(50, backend)
                arg = (info_level = -1, timesteps = :none)
                # Test basic version
                ws, states = simulate_reservoir(case; mode = :default, arg...)
                # Basic version (multiproc faked)
                @testset "PArray native" begin
                    for order in [:default]#, :symrcm]
                        for np in num_procs_to_test
                            ws_d, states_d = simulate_reservoir(case;
                                mode = :parray,
                                parray_arg = (np = np, order = order),
                                precond = :ilu0,
                                output_path = tempdir(),
                                arg...
                            )
                            compare_ws(ws, ws_d)
                            compare_states(states, states_d)
                        end
                    end
                end
                @testset "PArray MPI" begin
                    # MPI version
                    ws_m, states_m = simulate_reservoir(case; mode = :mpi, arg...)
                    compare_ws(ws, ws_m)
                    compare_states(states, states_m)
                end
            end
        end
    end
end
