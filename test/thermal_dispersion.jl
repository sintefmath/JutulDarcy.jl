using Test, Jutul, JutulDarcy
using JutulDarcy.StaticArrays: SVector

function thermal_dispersion_test_domain(; disp = nothing, hetero = true)
    dims = (8, 6, 1)
    if hetero
        dx = collect(range(20.0, 60.0, length = dims[1]))
        dy = collect(range(40.0, 15.0, length = dims[2]))
        g = CartesianMesh(dims, (dx, dy, [10.0]))
    else
        g = CartesianMesh(dims, (400.0, 300.0, 10.0))
    end
    nc = number_of_cells(g)
    poro = hetero ? collect(range(0.1, 0.4, length = nc)) : 0.25
    λ_f = hetero ? collect(range(0.4, 0.8, length = nc)) : 0.6
    domain = reservoir_domain(g,
        permeability = 0.5*si_unit(:darcy),
        porosity = poro,
        rock_thermal_conductivity = 2.0,
        fluid_thermal_conductivity = λ_f
    )
    if !isnothing(disp)
        JutulDarcy.add_thermal_dispersion!(domain, disp)
    end
    return domain
end

function thermal_dispersion_test_case(domain, sys; setup = missing, nstep = 5, kwarg...)
    bar = si_unit(:bar)
    day = si_unit(:day)
    nx, ny, = domain.representation.dims
    P = setup_vertical_well(domain, 1, 1, name = :Producer)
    I = setup_vertical_well(domain, nx, ny, name = :Injector)
    if sys == :geothermal
        model = setup_reservoir_model(domain, :geothermal; wells = [I, P], kwarg...)
        nph = 1
    else
        model = setup_reservoir_model(domain, sys; wells = [I, P], thermal = true, kwarg...)
        nph = number_of_phases(sys)
    end
    if !ismissing(setup)
        setup(model)
    end
    parameters = setup_parameters(model)
    mix = nph == 1 ? [1.0] : [1.0, 0.0]
    s0 = nph == 1 ? [1.0] : [0.6, 0.4]
    state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = s0, Temperature = 360.0)
    rate = 0.3*sum(pore_volume(model, parameters))/(nstep*10day)
    controls = Dict(
        :Injector => InjectorControl(TotalRateTarget(rate), mix, density = 1000.0, temperature = 300.0),
        :Producer => ProducerControl(BottomHolePressureTarget(100*bar))
    )
    forces = setup_reservoir_forces(model, control = controls)
    return JutulCase(model, fill(10day, nstep), forces, state0 = state0, parameters = parameters)
end

thermal_dispersion_final_temperature(case) = simulate_reservoir(case, info_level = -1).states[end][:Temperature]

single_phase = SinglePhaseSystem(AqueousPhase(), reference_density = 1000.0)
two_phase = ImmiscibleSystem((AqueousPhase(), LiquidPhase()), reference_densities = [1000.0, 800.0])

@testset "Dynamic fluid thermal conductivity" begin
    @testset "Reference conductivity reproduces static" begin
        for sys in (single_phase, two_phase)
            domain = thermal_dispersion_test_domain()
            T_static = thermal_dispersion_final_temperature(thermal_dispersion_test_case(domain, sys))
            case = thermal_dispersion_test_case(domain, sys,
                setup = m -> JutulDarcy.set_dynamic_fluid_thermal_conductivity!(m))
            @test haskey(case.parameters[:Reservoir], :FluidThermalHalfFaceWeights)
            T_dynamic = thermal_dispersion_final_temperature(case)
            @test T_dynamic ≈ T_static atol = 1e-8
        end
    end
    @testset "Weights give two-point transmissibility of ϕλ" begin
        domain = thermal_dispersion_test_domain()
        model = setup_reservoir_model(domain, single_phase, thermal = true)
        JutulDarcy.set_dynamic_fluid_thermal_conductivity!(model)
        prm = setup_parameters(model)[:Reservoir]
        w = prm[:FluidThermalHalfFaceWeights]
        θ_static = prm[:FluidThermalConductivities]
        λ = 2.0 .* prm[:FluidThermalConductivity] .* collect(range(0.5, 1.5, length = number_of_cells(domain)))'
        θ_ref = compute_face_trans(domain, domain[:porosity].*λ[1, :])
        N = domain[:neighbors]
        for face in axes(N, 2)
            l, r = N[1, face], N[2, face]
            f_l = JutulDarcy.fluid_thermal_conductivity_factor(λ, w, face, TPFA(l, r, 1), 1)
            f_r = JutulDarcy.fluid_thermal_conductivity_factor(λ, w, face, TPFA(r, l, -1), 1)
            @test f_l == f_r
            @test θ_static[1, face]*f_l ≈ θ_ref[face] rtol = 1e-12
        end
    end
    @testset "Derivatives of conduction and dispersion" begin
        θ_f = reshape([2.0], 1, 1)
        w = reshape([1.5, 3.0], 2, 1)
        θ_d = reshape([0.3, 0.7], 2, 1)
        λ(T) = 0.5 + 1e-3*T
        ρ(T) = 1000.0 - 0.3*T
        c(T) = 4000.0 + 2.0*T
        function conductive_flux(u, grad)
            row(f) = reshape(f.(u), 1, 2)
            state = (
                Temperature = u,
                FluidThermalConductivities = θ_f,
                FluidThermalHalfFaceWeights = w,
                FluidThermalConductivity = row(λ),
                ThermalDispersionTransmissibility = θ_d,
                PhaseMassDensities = row(ρ),
                FluidHeatCapacity = row(c)
            )
            θ = JutulDarcy.fluid_thermal_transmissibility(state, θ_f, 1, grad, 1) +
                JutulDarcy.thermal_dispersion_face_transmissibility(state, 1, grad, 1)
            return -θ*Jutul.gradient(u, grad)
        end
        function reference_flux(u)
            a_l, a_r = w[1]*λ(u[1]), w[2]*λ(u[2])
            d_l, d_r = θ_d[1]*ρ(u[1])*c(u[1]), θ_d[2]*ρ(u[2])*c(u[2])
            θ = θ_f[1]*a_l*a_r/(a_l + a_r) + d_l*d_r/(d_l + d_r)
            return -θ*(u[2] - u[1])
        end
        u = [330.0, 310.0]
        @test conductive_flux(u, TPFA(1, 2, 1)) ≈ reference_flux(u) rtol = 1e-12
        @test conductive_flux(u, TPFA(2, 1, -1)) ≈ -reference_flux(u) rtol = 1e-12
        derivatives = JutulDarcy.ForwardDiff.gradient(x -> conductive_flux(x, TPFA(1, 2, 1)), u)
        finite_difference = similar(derivatives)
        for i in eachindex(u)
            plus, minus = copy(u), copy(u)
            plus[i] += 1e-5
            minus[i] -= 1e-5
            finite_difference[i] = (reference_flux(plus) - reference_flux(minus))/2e-5
        end
        @test derivatives ≈ finite_difference rtol = 1e-7
        @test JutulDarcy.ForwardDiff.gradient(x -> conductive_flux(x, TPFA(2, 1, -1)), u) ≈ -derivatives rtol = 1e-12
    end
    @testset "Geothermal" begin
        domain = thermal_dispersion_test_domain(hetero = false)
        case = thermal_dispersion_test_case(domain, :geothermal, dynamic_conductivity = true)
        rmodel = reservoir_model(case.model)
        @test haskey(rmodel.secondary_variables, :FluidThermalConductivity)
        @test haskey(rmodel.parameters, :FluidThermalHalfFaceWeights)
        T = thermal_dispersion_final_temperature(case)
        @test all(isfinite, T)
    end
end

@testset "Thermal conductivity compatibility" begin
    domain = thermal_dispersion_test_domain(hetero = false)
    model = reservoir_model(setup_reservoir_model(domain, single_phase, thermal = true))
    prm = setup_parameters(model, RockThermalConductivities = 2.5, FluidThermalConductivities = 0.1)
    @test all(isequal(2.5), prm[:RockThermalConductivities])
    @test all(isequal(0.1), prm[:FluidThermalConductivities])
    @test !haskey(prm, :FluidThermalConductivity)
    domain[:fluid_thermal_conductivity] = fill(0.123, number_of_cells(domain))
    model = setup_reservoir_model(domain, :geothermal, update_reservoir = false)
    prm = setup_parameters(model)[:Reservoir]
    T_ref = compute_face_trans(domain, domain[:porosity].*domain[:fluid_thermal_conductivity])
    @test prm[:FluidThermalConductivities][1, :] ≈ T_ref
end

@testset "Thermal dispersion" begin
    disp = (10.0, 1.0)
    @testset "Dispersivity input" begin
        nc = 3
        @test JutulDarcy.thermal_dispersivity_values(2.0, nc) == fill(2.0, 2, nc)
        @test JutulDarcy.thermal_dispersivity_values((2.0, 1.0), nc) == repeat([2.0, 1.0], 1, nc)
        @test JutulDarcy.thermal_dispersivity_values([2.0, 1.0], nc) == repeat([2.0, 1.0], 1, nc)
        @test JutulDarcy.thermal_dispersivity_values([1.0, 2.0, 3.0], nc) == [1.0 2.0 3.0; 1.0 2.0 3.0]
        @test JutulDarcy.thermal_dispersivity_values([1.0 2.0 3.0; 0.1 0.2 0.3], nc) == [1.0 2.0 3.0; 0.1 0.2 0.3]
        @test_throws ArgumentError JutulDarcy.thermal_dispersivity_values([1.0, 2.0, 3.0, 4.0], nc)
        @test_throws ArgumentError JutulDarcy.thermal_dispersivity_values(-1.0, nc)
        domain = thermal_dispersion_test_domain(disp = disp)
        @test domain[:thermal_dispersivity] == repeat([10.0, 1.0], 1, number_of_cells(domain))
    end
    @testset "Half-face transmissibility" begin
        θ = JutulDarcy.dispersion_half_face_transmissibility
        α = reshape([10.0, 0.1], 2, 1)
        # Negative values on skewed geometry are set to zero
        @test θ(SVector(1.0, -2.0), SVector(1.0, 0.0), 1.0, SVector(1.0, 1.0), α, 1) == 0.0
        # The normal velocity comes from the face flux, even if the cell velocity vanishes
        @test θ(SVector(0.5, 0.0), SVector(2.0, 0.0), 3.0, SVector(0.0, 0.0), α, 1) ≈ 10.0*1.5*0.5*2.0/0.25
    end
    @testset "Zero dispersivity" begin
        T_none = thermal_dispersion_final_temperature(
            thermal_dispersion_test_case(thermal_dispersion_test_domain(), single_phase))
        T_zero = thermal_dispersion_final_temperature(
            thermal_dispersion_test_case(thermal_dispersion_test_domain(disp = 0.0), single_phase))
        @test T_zero == T_none
    end
    @testset "Geothermal with $(dyn ? "dynamic" : "static") conductivity" for dyn in (false, true)
        domain = thermal_dispersion_test_domain(disp = disp)
        case = thermal_dispersion_test_case(domain, :geothermal, dynamic_conductivity = dyn)
        rmodel = reservoir_model(case.model)
        @test rmodel.secondary_variables[:FluidHeatCapacity] isa JutulDarcy.PressureTemperatureDependentVariable
        T = thermal_dispersion_final_temperature(case)
        T_none = thermal_dispersion_final_temperature(
            thermal_dispersion_test_case(thermal_dispersion_test_domain(), :geothermal, dynamic_conductivity = dyn))
        @test all(isfinite, T)
        @test maximum(abs, T - T_none) > 0.1
    end
    @testset "Two-phase flux is conservative" begin
        domain = thermal_dispersion_test_domain(disp = disp)
        case = thermal_dispersion_test_case(domain, two_phase,
            setup = m -> JutulDarcy.set_dynamic_fluid_thermal_conductivity!(m))
        sim = Simulator(case)
        simulate!(sim, case.dt[1:2], forces = case.forces, info_level = -1)
        rmodel = reservoir_model(case.model)
        state = sim.storage[:Reservoir].state
        @test size(state.ThermalDispersionTransmissibility, 1) == 4
        @test maximum(state.ThermalDispersionTransmissibility) > 0
        N = domain[:neighbors]
        for face in axes(N, 2)
            l, r = N[1, face], N[2, face]
            q_l = JutulDarcy.thermal_heat_flux(face, state, rmodel, TPFA(l, r, 1), SPU(l, r), Jutul.DefaultFlux())
            q_r = JutulDarcy.thermal_heat_flux(face, state, rmodel, TPFA(r, l, -1), SPU(r, l), Jutul.DefaultFlux())
            @test value(q_l) ≈ -value(q_r) rtol = 1e-12
        end
    end
    @testset "1D dispersion matches equivalent conductivity" begin
        bar = si_unit(:bar)
        day = si_unit(:day)
        nx = 30
        A = 100.0
        q = 1e-4
        α_L = 5.0
        function run_1d(α, extra_conductivity)
            g = CartesianMesh((nx, 1, 1), (300.0, 10.0, 10.0))
            ϕ = 0.2
            domain = reservoir_domain(g,
                permeability = 1.0*si_unit(:darcy),
                porosity = ϕ,
                rock_thermal_conductivity = 2.0 + extra_conductivity/(1.0 - ϕ),
                fluid_thermal_conductivity = 0.6
            )
            if α > 0
                JutulDarcy.add_thermal_dispersion!(domain, (α, 0.0))
            end
            I = setup_well(domain, [(1, 1, 1)], name = :I, simple_well = true)
            P = setup_well(domain, [(nx, 1, 1)], name = :P, simple_well = true)
            model = setup_reservoir_model(domain, single_phase, wells = [I, P], thermal = true)
            state0 = setup_reservoir_state(model, Pressure = 100bar, Temperature = 360.0)
            controls = Dict(
                :I => InjectorControl(TotalRateTarget(q), [1.0], density = 1000.0, temperature = 300.0),
                :P => ProducerControl(BottomHolePressureTarget(90bar))
            )
            forces = setup_reservoir_forces(model, control = controls)
            case = JutulCase(model, fill(5day, 20), forces, state0 = state0)
            return thermal_dispersion_final_temperature(case)
        end
        ρc = 1000.0*JutulDarcy.default_value(nothing, JutulDarcy.ComponentHeatCapacity())
        T_dispersion = run_1d(α_L, 0.0)
        T_conductivity = run_1d(0.0, ρc*α_L*q/A)
        @test T_dispersion ≈ T_conductivity atol = 0.05
    end
end
