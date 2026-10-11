using Jutul, JutulDarcy, JPhreeqc, Test
using LinearAlgebra: I
import JutulDarcy.Tracers as Tracers
const PhreeqcFD = JutulDarcy.ForwardDiff

const PHREEQC_VALIDATION_INPUT = """
SOLUTION 1
    temp 25
    pH 9
    units mol/kgw
    Na 0.01
    Cl 0.01
    Ca 0.002
    C(4) 0.004
SOLUTION 0
    temp 25
    pH 7
    units mol/kgw
    Na 0.01
    Cl 0.01
    Ca 0.0001
    C(4) 0.0002
END
"""

function phreeqc_test_model(nc = 4; phases = ["Calcite"], batch_size = 1, kwarg...)
    mesh = CartesianMesh((nc, 1), (Float64(nc), 1.0))
    domain = reservoir_domain(mesh, permeability = 1e-12, porosity = 0.2)
    model = setup_reservoir_model(domain, SinglePhaseSystem(AqueousPhase()))
    rmodel = reservoir_model(model)
    rmodel.parameters[:Temperature] = JutulDarcy.Temperature()
    setup_phreeqc_model!(model, PHREEQC_VALIDATION_INPUT;
        equilibrium_phases = phases, batch_size = batch_size, kwarg...)
    return model
end

function phreeqc_test_native_threads()
    id = RM_Create(8, 4)
    id >= 0 || error("Could not probe native thread support.")
    try
        return RM_GetThreadCount(id)
    finally
        RM_Destroy(id)
    end
end

@testset "PHREEQC task-local prepared workers" begin
    model = phreeqc_test_model(4)
    try
        ext = Base.get_extension(JutulDarcy, :JutulDarcyJPhreeqcExt)
        chemistry = reservoir_model(model).secondary_variables[:PhreeqcPhaseSplit]
        @test !haskey(task_local_storage(), chemistry.task_key)
        ext.with_chemistry_worker(chemistry) do worker
            @test worker in chemistry.workers
            @test ext.chemistry_worker(chemistry) === worker
            ext.with_chemistry_worker(chemistry) do nested
                @test nested === worker
                yield()
                @test ext.chemistry_worker(chemistry) === worker
            end
        end
        @test !haskey(task_local_storage(), chemistry.task_key)
        @test_throws ErrorException ext.with_chemistry_worker(chemistry) do worker
            error("Test worker lease cleanup.")
        end
        @test !haskey(task_local_storage(), chemistry.task_key)
        @test Base.n_avail(chemistry.available) == length(chemistry.workers)
    finally
        close_phreeqc!(model)
    end
end

@testset "PHREEQC heterogeneous batches and padding" begin
    batched = phreeqc_test_model(5, batch_size = 3)
    serial = phreeqc_test_model(5, threads = 1)
    try
        c = phreeqc_initial_concentrations(batched)
        n = length(c)
        totals = repeat(c, 1, 5)
        ca = findfirst(==("Ca"), phreeqc_component_names(batched))
        totals[ca, :] .*= range(1.0, 1.2, length = 5)
        function evaluate(model, flat)
            rmodel = reservoir_model(model)
            state = (TracerConcentrations = reshape(flat, n, 5),
                Pressure = collect(range(101325.0, 2e5, length = 5)),
                Temperature = collect(range(298.15, 323.15, length = 5)),
                Saturations = ones(1, 5),
                PhaseMassDensities = reshape(collect(range(1000.0, 1010.0, length = 5)), 1, 5))
            chemistry = rmodel.secondary_variables[:PhreeqcPhaseSplit]
            out = zeros(eltype(flat), Jutul.values_per_entity(rmodel, chemistry), 5)
            Jutul.update_secondary_variable!(out, chemistry, rmodel, state, 1:5)
            return vec(out[1:n, :])
        end
        flat = vec(totals)
        @test evaluate(batched, flat) ≈ evaluate(serial, flat) atol = 1e-12
        jb = PhreeqcFD.jacobian(x -> evaluate(batched, x), flat)
        js = PhreeqcFD.jacobian(x -> evaluate(serial, x), flat)
        @test jb ≈ js rtol = 3e-3 atol = 3e-5
        @test all(iszero, jb[1:n, (n + 1):end])
    finally
        close_phreeqc!(batched)
        close_phreeqc!(serial)
    end
end

function phreeqc_test_split(model, x)
    rmodel = reservoir_model(model)
    chemistry = rmodel.secondary_variables[:PhreeqcPhaseSplit]
    n = length(chemistry.components)
    nc = number_of_cells(rmodel.domain)
    state = (TracerConcentrations = repeat(reshape(x[1:n], n, 1), 1, nc),
        Pressure = fill(x[n + 1], nc), Temperature = fill(x[n + 2], nc),
        Saturations = ones(eltype(x), 1, nc),
        PhaseMassDensities = reshape(fill(x[n + 3], nc), 1, nc))
    out = zeros(eltype(x), Jutul.values_per_entity(rmodel, chemistry), nc)
    Jutul.update_secondary_variable!(out, chemistry, rmodel, state, 1:nc)
    return out[1:n, :]
end

function independent_calcite_reference()
    id = IPQ_CreateIPhreeqc()
    try
        @test IPQ_LoadDatabase(id, database_path("phreeqc.dat")) == 0
        input = PHREEQC_VALIDATION_INPUT * """
        SELECTED_OUTPUT 1
            -reset false
            -high_precision true
        USER_PUNCH 1
            -headings Ca C Calcite
            -start
            10 PUNCH TOTMOLE("Ca"), TOTMOLE("C"), EQUI("Calcite")
            -end
        EQUILIBRIUM_PHASES 1
            Calcite 0 0
        USE solution 1
        END
        """
        @test IPQ_RunString(id, input) == 0
        _, values = getIPhreeqcSelectedOutput(id, 1)
        return values[end, :]
    finally
        IPQ_DestroyIPhreeqc(id)
    end
end

@testset "PHREEQC conservative equilibrium split" begin
    model = phreeqc_test_model()
    serial = phreeqc_test_model(4, threads = 1, reuse_speciation = false)
    try
        names = phreeqc_component_names(model)
        c = phreeqc_initial_concentrations(model)
        x = [c; 101325.0; 298.15; 1000.0]
        y = phreeqc_test_split(model, x)
        @test y ≈ phreeqc_test_split(serial, x) rtol = 1e-12 atol = 1e-12
        @test all(isfinite, y)
        ca = findfirst(==("Ca"), names)
        carbon = findfirst(==("C"), names)
        oxygen = findfirst(==("O"), names)
        hydrogen = findfirst(==("H"), names)
        molar_mass = phreeqc_component_molar_masses(model)
        ref_ca, ref_c, ref_mineral = independent_calcite_reference()
        carrier_mass = 0.002*molar_mass[ca]/c[ca]
        @test y[ca, 1] ≈ ref_ca*molar_mass[ca]/carrier_mass rtol = 2e-6
        @test y[carbon, 1] ≈ ref_c*molar_mass[carbon]/carrier_mass rtol = 2e-6
        mineral_ca = c[ca] - y[ca, 1]
        mineral_c = c[carbon] - y[carbon, 1]
        @test mineral_ca > 0
        @test mineral_ca ≈ ref_mineral*molar_mass[ca]/carrier_mass rtol = 2e-6
        @test mineral_ca/molar_mass[ca] ≈ mineral_c/molar_mass[carbon] rtol = 2e-6
        @test (c[oxygen] - y[oxygen, 1])/molar_mass[oxygen] ≈ 3*mineral_c/molar_mass[carbon] rtol = 2e-6
        @test abs(c[hydrogen] - y[hydrogen, 1]) < 1e-11
        # A different evaluation must not advance mineral history. Return to
        # the original state and demand the original answer after a cache miss.
        altered = copy(x)
        altered[ca] *= 1.1
        phreeqc_test_split(model, altered)
        restored = phreeqc_test_split(model, x)
        @test restored ≈ y rtol = 1e-12 atol = 1e-12
        @test phreeqc_test_split(model, x) == restored
        @test deepcopy(reservoir_model(model).secondary_variables[:PhreeqcPhaseSplit]) ===
            reservoir_model(model).secondary_variables[:PhreeqcPhaseSplit]
        # The well concentration hook must return the mobile fraction.
        tracer = reservoir_model(model).equations[:tracers].flux_type.tracers[ca]
        state = (TracerConcentrations = repeat(c, 1, 4), PhreeqcAqueousConcentrations = y)
        @test Tracers.tracer_phase_concentration(tracer, reservoir_model(model), state, ca, 1, 1) == y[ca, 1]
        close_phreeqc!(model)
        close_phreeqc!(model)
        @test_throws ErrorException phreeqc_test_split(model, x)
    finally
        close_phreeqc!(model)
        close_phreeqc!(serial)
    end
end

@testset "PHREEQC phase-split partial derivatives" begin
    model = phreeqc_test_model(2, difference = :central, relative_step = 2e-5)
    try
        c = phreeqc_initial_concentrations(model)
        n = length(c)
        x = [c; 2e5; 310.0; 1010.0]
        f = x -> phreeqc_test_split(model, x)[:, 1]
        jacobian = PhreeqcFD.jacobian(f, x)
        @test all(isfinite, jacobian)
        for j in eachindex(x)
            scale = j <= n ? max(abs(x[j]), 1e-6) : (j == n + 1 ? 1e5 : 1000.0)
            h = 1e-4*scale
            plus = copy(x); plus[j] += h
            minus = copy(x); minus[j] -= h
            ref = (f(plus) - f(minus))/(2*h)
            @test jacobian[:, j] ≈ ref rtol = 3e-3 atol = 3e-6
        end
        # Every conserved inventory splits into mobile + immobile, including
        # derivatives at fixed pressure, temperature and density.
        joint = x -> begin
            aqueous = f(x)
            return [aqueous; x[1:n] - aqueous]
        end
        split_jacobian = PhreeqcFD.jacobian(joint, x)
        identity = zeros(n, n + 3)
        for i in 1:n
            identity[i, i] = 1
        end
        @test split_jacobian[1:n, :] + split_jacobian[(n + 1):end, :] ≈ identity atol = 1e-12
        ca = findfirst(==("Ca"), phreeqc_component_names(model))
        carbon = findfirst(==("C"), phreeqc_component_names(model))
        @test abs(jacobian[ca, carbon]) > 1e-3
        @test abs(jacobian[ca, n + 2]) > 1e-10
    finally
        close_phreeqc!(model)
    end
end

@testset "PHREEQC mineral volume closure and derivatives" begin
    model = phreeqc_test_model(1, difference = :central, relative_step = 2e-5,
        mineral_molar_volumes = [3.69e-5], reference_mineral_volume = 1e-3)
    try
        rmodel = reservoir_model(model)
        chemistry = rmodel.secondary_variables[:PhreeqcPhaseSplit]
        volume = rmodel.secondary_variables[:FluidVolume]
        c = phreeqc_initial_concentrations(model)
        names = phreeqc_component_names(model)
        masses = phreeqc_component_molar_masses(model)
        n = length(c)
        for (name, atoms) in (("Ca", 1), ("C", 1), ("O", 3))
            i = findfirst(==(name), names)
            c[i] += 0.2*atoms*masses[i]
        end
        function evaluate(x)
            state = (TracerConcentrations = reshape(x[1:n], n, 1),
                Pressure = [x[n + 1]], Temperature = [x[n + 2]],
                PhaseMassDensities = reshape([x[n + 3]], 1, 1),
                Saturations = reshape([x[n + 4]], 1, 1))
            split = zeros(eltype(x), n + 2, 1)
            Jutul.update_secondary_variable!(split, chemistry, rmodel, state, 1:1)
            pvstate = (; state..., StaticFluidVolume = [x[n + 5]],
                PhreeqcReferenceMineralVolume = [x[n + 6]], PhreeqcMineralVolume = [split[end, 1]])
            pv = zeros(eltype(x), 1)
            Jutul.update_secondary_variable!(pv, volume, rmodel, pvstate, 1:1)
            return [vec(split); pv]
        end
        x = [c; 2e5; 310.0; 1010.0; 0.7; 0.2; 1e-3]
        y = evaluate(x)
        @test y[n + 1] > 0.19 # mol calcite/kg carrier
        @test y[n + 2] ≈ 3.69e-5*y[n + 1] rtol = 1e-12
        solid_volume = x[n + 3]*x[n + 4]*y[end]*y[n + 2]
        @test y[end] + solid_volume ≈ x[n + 5] + x[n + 6] rtol = 1e-12
        @test y[end] < x[n + 5]
        # Complete chain rule through phase amounts, density, saturation,
        # reference pore volume and pre-existing mineral volume.
        jacobian = PhreeqcFD.jacobian(evaluate, x)
        for j in eachindex(x)
            scale = j <= n ? 1e-6 : (j == n + 1 ? 1e5 : (j == n + 2 ? 300.0 : 1.0))
            h = 1e-4*max(abs(x[j]), scale)
            plus = copy(x); plus[j] += h
            minus = copy(x); minus[j] -= h
            ref = (evaluate(plus) - evaluate(minus))/(2*h)
            @test jacobian[:, j] ≈ ref rtol = 4e-3 atol = 3e-6
        end
        @test all(iszero, jacobian[1:(n + 2), n + 3]) # exact density independence
        @test jacobian[end, n + 3] < 0
        @test jacobian[end, n + 4] < 0
        # Dissolving the initial solid opens pore space and releases the
        # specified reference solid volume, with no timestep history.
        dilute = copy(x)
        dilute[1:n] = phreeqc_initial_concentrations(model, solution = 0)
        dissolved = evaluate(dilute)
        @test dissolved[n + 1] == 0
        @test dissolved[end] ≈ x[n + 5] + x[n + 6]
        @test dissolved[end] > y[end]
        @test evaluate(x) ≈ y rtol = 1e-12 atol = 1e-12
    finally
        close_phreeqc!(model)
    end
end

@testset "PHREEQC exact inactive and spectator derivatives" begin
    plain = phreeqc_test_model(1, phases = String[])
    mineral = phreeqc_test_model(1)
    try
        for model in (plain, mineral)
            c = phreeqc_initial_concentrations(model, solution = 0)
            x = [c; 101325.0; 298.15; 1000.0]
            f = x -> phreeqc_test_split(model, x)[:, 1]
            jacobian = PhreeqcFD.jacobian(f, x)
            @test f(x) == reshape(c, length(c), 1)[:, 1]
            @test jacobian[:, 1:length(c)] == Matrix{Float64}(I, length(c), length(c))
            @test all(iszero, jacobian[:, (length(c) + 1):end])
        end
        c = phreeqc_initial_concentrations(mineral)
        n = length(c)
        jacobian = PhreeqcFD.jacobian(x -> phreeqc_test_split(mineral, x)[:, 1],
            [c; 101325.0; 298.15; 1000.0])
        names = phreeqc_component_names(mineral)
        for name in ("H", "Na", "Cl", "Charge")
            i = findfirst(==(name), names)
            @test jacobian[i, :] == Float64.(collect(1:(n + 3)) .== i)
        end
    finally
        close_phreeqc!(plain)
        close_phreeqc!(mineral)
    end
end

@testset "PHREEQC reactive and compressible pore volume" begin
    mesh = CartesianMesh((1, 1), (1.0, 1.0))
    model = setup_reservoir_model(reservoir_domain(mesh, porosity = 0.2), SinglePhaseSystem(AqueousPhase()))
    rmodel = reservoir_model(model)
    JutulDarcy.set_rock_compressibility!(model, compressibility = 1e-8, reference_pressure = 101325.0)
    setup_phreeqc_model!(model, PHREEQC_VALIDATION_INPUT, equilibrium_phases = ["Calcite"])
    fixed = phreeqc_test_model(1, update_pore_volume = false)
    try
        c = phreeqc_initial_concentrations(model)
        state0 = setup_reservoir_state(model, Pressure = 2e5, TracerConcentrations = c)
        state = Jutul.evaluate_all_secondary_variables(rmodel, state0[:Reservoir])
        @test state[:PhreeqcUnreactedFluidVolume][1] > 0.2
        @test state[:FluidVolume][1] < state[:PhreeqcUnreactedFluidVolume][1]
        occupied = state[:PhreeqcMineralVolume][1]*state[:PhaseMassDensities][1]*state[:FluidVolume][1]
        @test state[:FluidVolume][1] + occupied ≈ state[:PhreeqcUnreactedFluidVolume][1] rtol = 1e-12
        cfixed = phreeqc_initial_concentrations(fixed)
        state0 = setup_reservoir_state(fixed, Pressure = 101325.0, TracerConcentrations = cfixed)
        state = Jutul.evaluate_all_secondary_variables(reservoir_model(fixed), state0[:Reservoir])
        @test state[:PhreeqcMineralVolume][1] > 0
        @test state[:FluidVolume][1] == 0.2
    finally
        close_phreeqc!(model)
        close_phreeqc!(fixed)
    end
end

@testset "PHREEQC distinct mineral amounts and database volumes" begin
    model = phreeqc_test_model(1, phases = ["Calcite", "Aragonite"])
    try
        c = phreeqc_initial_concentrations(model)
        rmodel = reservoir_model(model)
        state0 = setup_reservoir_state(model, Pressure = 101325.0, TracerConcentrations = c)
        state = Jutul.evaluate_all_secondary_variables(rmodel, state0[:Reservoir])
        mineral = state[:PhreeqcMineralConcentrations][:, 1]
        @test mineral[1] > 0
        @test mineral[2] == 0 # Polymorphs share a formula but remain distinct phases.
        @test state[:PhreeqcMineralVolume][1] ≈ mineral[1]*3.69e-5 rtol = 1e-12
        chemistry = rmodel.secondary_variables[:PhreeqcPhaseSplit]
        @test first(chemistry.workers).molar_volumes ≈ [3.69e-5, 3.404e-5] rtol = 1e-12
    finally
        close_phreeqc!(model)
    end
    mesh = CartesianMesh((1, 1), (1.0, 1.0))
    missing_volume = """
    PHASES
    NoVolumeCalcite
        CaCO3 = Ca+2 + CO3-2
        log_k -8.48
    """
    model = setup_reservoir_model(reservoir_domain(mesh), SinglePhaseSystem(AqueousPhase()))
    @test_throws ArgumentError setup_phreeqc_model!(model,
        missing_volume * PHREEQC_VALIDATION_INPUT, equilibrium_phases = ["NoVolumeCalcite"])
    # A failed setup must leave the flow model usable for another setup.
    setup_phreeqc_model!(model, missing_volume * PHREEQC_VALIDATION_INPUT,
        equilibrium_phases = ["NoVolumeCalcite"], mineral_molar_volumes = [3.69e-5])
    try
        @test length(phreeqc_initial_concentrations(model)) > 0
    finally
        close_phreeqc!(model)
    end
end

@testset "PHREEQC transport conserves total inventory" begin
    for threading in (:julia, :native)
        threads = threading == :native ? phreeqc_test_native_threads() : Threads.nthreads()
        batch_size = threading == :native ? nothing : 1
        model = phreeqc_test_model(4, threading = threading, threads = threads, batch_size = batch_size)
        try
            c = phreeqc_initial_concentrations(model)
            totals = repeat(c, 1, 4)
            ca = findfirst(==("Ca"), phreeqc_component_names(model))
            totals[ca, 1] *= 1.2
            state0 = setup_reservoir_state(model,
                Pressure = [2e5, 1.8e5, 1.5e5, 1.1e5], TracerConcentrations = totals)
            rmodel = reservoir_model(model)
            initial = Jutul.evaluate_all_secondary_variables(rmodel, state0[:Reservoir])
            sim = Simulator(rmodel, state0 = state0[:Reservoir], parameters = setup_parameters(rmodel))
            cfg = simulator_config(sim, info_level = -1)
            states, reports = simulate(sim, [0.01, 0.01], config = cfg)
            @test length(states) == 2
            before = vec(sum(initial[:TracerMasses], dims = 2))
            after = vec(sum(states[end][:TracerMasses], dims = 2))
            @test after ≈ before rtol = 2e-7 atol = 1e-9
            @test sum(states[end][:TotalMasses]) ≈ sum(initial[:TotalMasses]) rtol = 2e-7
            @test maximum(abs.(states[end][:TracerConcentrations] .- totals)) > 1e-10
        finally
            close_phreeqc!(model)
        end
    end
end

@testset "PHREEQC in Jutul's threaded secondary loop" begin
    # Exceed ParallelCSRContext's default minimum batch size, so Jutul
    # partitions all secondary variables before entering the chemistry update.
    model = phreeqc_test_model(2048, batch_size = 32)
    try
        rmodel = reservoir_model(model)
        c = phreeqc_initial_concentrations(model)
        state0 = setup_reservoir_state(model, Pressure = 101325.0, TracerConcentrations = c)
        state = Jutul.evaluate_all_secondary_variables(rmodel, state0[:Reservoir])
        @test state[:PhreeqcAqueousConcentrations][:, 1] ≈ state[:PhreeqcAqueousConcentrations][:, end] atol = 1e-12
        @test all(isfinite, state[:PhreeqcPhaseSplit])
        occupied = state[:PhreeqcMineralVolume] .* state[:PhaseMassDensities][1, :] .* state[:FluidVolume]
        @test state[:FluidVolume] + occupied ≈ fill(0.2, 2048) rtol = 1e-12
    finally
        close_phreeqc!(model)
    end
end

@testset "PHREEQC validation and native error recovery" begin
    model = phreeqc_test_model(1)
    try
        x = [phreeqc_initial_concentrations(model); 101325.0; 298.15; 1000.0]
        ca = findfirst(==("Ca"), phreeqc_component_names(model))
        invalid = copy(x); invalid[ca] = -1.0
        @test_throws ErrorException phreeqc_test_split(model, invalid)
        @test all(isfinite, phreeqc_test_split(model, x))
        @test_throws ArgumentError setup_phreeqc_model!(model, PHREEQC_VALIDATION_INPUT)
        @test_throws ArgumentError phreeqc_initial_concentrations(model, solution = 999)
    finally
        close_phreeqc!(model)
    end
    mesh = CartesianMesh((1, 1), (1.0, 1.0))
    plain = setup_reservoir_model(reservoir_domain(mesh), SinglePhaseSystem(AqueousPhase()))
    @test_throws ArgumentError setup_phreeqc_model!(plain, "SOLUTION 1\nEND\nKINETICS 1\nEND")
    @test_throws ArgumentError setup_phreeqc_model!(plain, PHREEQC_VALIDATION_INPUT, batch_size = 0)
    @test_throws ArgumentError setup_phreeqc_model!(plain, PHREEQC_VALIDATION_INPUT, threading = :invalid)
    @test_throws ArgumentError setup_phreeqc_model!(plain, PHREEQC_VALIDATION_INPUT, rebalance_fraction = 1.1)
    @test_throws ErrorException setup_phreeqc_model!(plain, "SOLUTION 1\nNa invalid_concentration\nEND\n")
end

@testset "PHREEQC native cell threading and pore-volume derivatives" begin
    available_threads = phreeqc_test_native_threads()
    serial = phreeqc_test_model(8, threads = 1, batch_size = 3,
        difference = :central, reference_mineral_volume = 1e-3)
    try
        for threads in unique([1, available_threads])
            model = phreeqc_test_model(8, threading = :native, threads = threads,
                batch_size = nothing, difference = :central, reference_mineral_volume = 1e-3)
            try
                rmodel = reservoir_model(model)
                chemistry = rmodel.secondary_variables[:PhreeqcPhaseSplit]
                @test length(chemistry.workers) == 1
                @test chemistry.batch_size == 8
                @test RM_GetThreadCount(first(chemistry.workers).id) == threads
                c = phreeqc_initial_concentrations(model)
                ca = findfirst(==("Ca"), phreeqc_component_names(model))
                n = length(c)
                function evaluate(model, x; ix = 1:8)
                    rmodel = reservoir_model(model)
                    totals = repeat(x[1:n], 1, 8)
                    totals[ca, :] .*= range(1.0, 1.2, length = 8)
                    state = (TracerConcentrations = totals, Pressure = fill(x[n + 1], 8),
                        Temperature = x[n + 2] .+ collect(range(0.0, 5.0, length = 8)),
                        PhaseMassDensities = fill(x[n + 3], 1, 8), Saturations = fill(x[n + 4], 1, 8))
                    split = zeros(eltype(x), n + 2, 8)
                    Jutul.update_secondary_variable!(split, rmodel.secondary_variables[:PhreeqcPhaseSplit],
                        rmodel, state, ix)
                    pvstate = (; state..., StaticFluidVolume = fill(x[n + 5], 8),
                        PhreeqcReferenceMineralVolume = fill(x[n + 6], 8),
                        PhreeqcMineralVolume = vec(split[end, :]))
                    pv = zeros(eltype(x), 8)
                    Jutul.update_secondary_variable!(pv, rmodel.secondary_variables[:FluidVolume], rmodel, pvstate, ix)
                    return vcat(split, transpose(pv))
                end
                x = [c; 2e5; 310.0; 1010.0; 0.7; 0.2; 1e-3]
                value = evaluate(model, x)
                @test value ≈ evaluate(serial, x) rtol = 1e-9 atol = 1e-12
                jacobian = PhreeqcFD.jacobian(x -> vec(evaluate(model, x)), x)
                reference = PhreeqcFD.jacobian(x -> vec(evaluate(serial, x)), x)
                @test jacobian ≈ reference rtol = 3e-3 atol = 3e-6
                for j in eachindex(x)
                    scale = j <= n ? 1e-6 : (j == n + 1 ? 1e5 : (j == n + 2 ? 300.0 : 1.0))
                    h = 1e-4*max(abs(x[j]), scale)
                    plus = copy(x); plus[j] += h
                    minus = copy(x); minus[j] -= h
                    independent = vec(evaluate(model, plus) - evaluate(model, minus))/(2*h)
                    @test jacobian[:, j] ≈ independent rtol = 4e-3 atol = 3e-6
                end
                @test value[end, :] + x[n + 3]*x[n + 4]*value[end, :].*value[n + 2, :] ≈
                    fill(x[n + 5] + x[n + 6], 8) rtol = 1e-12
                # Reverse and overlapping partial ranges retain cell ordering
                # and never share the single OpenMP module concurrently.
                chemistry.cache_valid .= false
                @test evaluate(model, x, ix = 8:-1:1) ≈ value rtol = 1e-9 atol = 1e-12
                chemistry.cache_valid .= false
                a, b = similar(value), similar(value)
                @sync begin
                    Threads.@spawn copyto!(a, evaluate(model, x, ix = 1:4))
                    Threads.@spawn copyto!(b, evaluate(model, x, ix = 5:8))
                end
                @test a[:, 1:4] ≈ value[:, 1:4] rtol = 1e-9 atol = 1e-12
                @test b[:, 5:8] ≈ value[:, 5:8] rtol = 1e-9 atol = 1e-12
                chemistry.cache_valid .= false
                @sync begin
                    Threads.@spawn copyto!(a, evaluate(model, x))
                    Threads.@spawn copyto!(b, evaluate(model, x))
                end
                @test a ≈ value rtol = 1e-9 atol = 1e-12
                @test b ≈ value rtol = 1e-9 atol = 1e-12
                @test Base.n_avail(chemistry.available) == 1
            finally
                close_phreeqc!(model)
            end
        end
    finally
        close_phreeqc!(serial)
    end
end

@testset "PHREEQC partial ranges and shared runtime" begin
    model = phreeqc_test_model(8)
    try
        c = phreeqc_initial_concentrations(model)
        x = [c; 101325.0; 298.15; 1000.0]
        expected = phreeqc_test_split(model, x)
        rmodel = reservoir_model(model)
        chemistry = rmodel.secondary_variables[:PhreeqcPhaseSplit]
        state = (TracerConcentrations = repeat(c, 1, 8), Pressure = fill(101325.0, 8),
            Saturations = ones(1, 8),
            Temperature = fill(298.15, 8), PhaseMassDensities = fill(1000.0, 1, 8))
        chemistry.cache_valid .= false
        out = zeros(Jutul.values_per_entity(rmodel, chemistry), 8)
        @sync for ix in (1:4, 5:8)
            Threads.@spawn Jutul.update_secondary_variable!(out, chemistry, rmodel, state, ix)
        end
        @test out[1:length(c), :] ≈ expected atol = 1e-12
        # Overlapping calls on a shared copied-model runtime serialize only
        # the affected batches and return independently correct outputs.
        a = similar(out)
        b = similar(a)
        chemistry.cache_valid .= false
        @sync begin
            Threads.@spawn Jutul.update_secondary_variable!(a, chemistry, rmodel, state, 1:8)
            Threads.@spawn Jutul.update_secondary_variable!(b, chemistry, rmodel, state, 1:8)
        end
        @test a[1:length(c), :] ≈ expected atol = 1e-12
        @test b[1:length(c), :] ≈ expected atol = 1e-12
    finally
        close_phreeqc!(model)
    end
end

@testset "PHREEQC aqueous carrier in multiphase flow" begin
    mesh = CartesianMesh((2, 1), (2.0, 1.0))
    system = ImmiscibleSystem((LiquidPhase(), AqueousPhase()), reference_densities = [800.0, 1000.0])
    model = setup_reservoir_model(reservoir_domain(mesh, porosity = 0.2), system)
    setup_phreeqc_model!(model, PHREEQC_VALIDATION_INPUT, equilibrium_phases = ["Calcite"])
    try
        c = phreeqc_initial_concentrations(model)
        ca = findfirst(==("Ca"), phreeqc_component_names(model))
        totals = repeat(c, 1, 2)
        totals[ca, 1] *= 1.1
        state0 = setup_reservoir_state(model, Pressure = 101325.0,
            Saturations = [0.3, 0.7], TracerConcentrations = totals)
        rmodel = reservoir_model(model)
        evaluated = Jutul.evaluate_all_secondary_variables(rmodel, state0[:Reservoir])
        state = Jutul.JutulStorage(evaluated)
        carrier = state.FluidVolume .* state.Saturations[2, :] .* state.PhaseMassDensities[2, :]
        @test state.TracerMasses ≈ totals .* reshape(carrier, 1, 2)
        tracer = rmodel.equations[:tracers].flux_type.tracers[ca]
        upwind = Jutul.SPU(1, 2)
        @test Tracers.tracer_total_mass_flux(tracer, rmodel, state, [100.0, 2.0], ca, upwind) ≈
            2*state.PhreeqcAqueousConcentrations[ca, 1]
        @test Tracers.tracer_total_mass_flux(tracer, rmodel, state, [100.0, -2.0], ca, upwind) ≈
            -2*state.PhreeqcAqueousConcentrations[ca, 2]
        state.Saturations[2, 1] = 0.0
        chemistry = rmodel.secondary_variables[:PhreeqcPhaseSplit]
        @test_throws ErrorException Jutul.update_secondary_variable!(state.PhreeqcPhaseSplit,
            chemistry, rmodel, state, 1:2)
    finally
        close_phreeqc!(model)
    end
end
