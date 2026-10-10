using Jutul, JutulDarcy, MultiComponentFlash, Test

function diffusion_test_model(domain, system)
    model = setup_reservoir_model(domain, system)
    parameters = setup_parameters(model)
    return reservoir_model(model), parameters[:Reservoir]
end

@testset "Molar diffusion values and derivatives" begin
    mw = [18.015, 16.043, 2.016]./1000
    S = [0.2 0.4; 0.8 0.6]
    rho = [1000.0 1100.0; 5.0 7.0]
    Dl = reshape([1e-9, 2e-9, 4e-9], 3, 1)
    Dv = 2Dl
    grad = TPFA(1, 2, 1)
    # Independent mole fractions in both phases/cells, with the last fraction
    # recovered from the unit-sum constraint.
    u = [0.7, 0.2, 0.6, 0.1, 0.1, 0.6, 0.2, 0.5]
    function fractions(u)
        X = [u[1] u[3]; u[2] u[4]; 1-u[1]-u[2] 1-u[3]-u[4]]
        Y = [u[5] u[7]; u[6] u[8]; 1-u[5]-u[6] 1-u[7]-u[8]]
        return X, Y
    end
    function flux(u, grad = grad, saturations = S)
        X, Y = fractions(u)
        liquid = mw.*X./sum(mw.*X; dims = 1)
        vapor = mw.*Y./sum(mw.*Y; dims = 1)
        q = JutulDarcy.SVector{3}(zero(eltype(u)), zero(eltype(u)), zero(eltype(u)))
        q = JutulDarcy.add_phase_diffusive_component_flux(q, Dl, 1, grad,
            saturations, rho, liquid, 1, mw, Val(3))
        return JutulDarcy.add_phase_diffusive_component_flux(q, Dv, 1, grad,
            saturations, rho, vapor, 2, mw, Val(3))
    end
    function reference_flux(u)
        X, Y = fractions(u)
        result = zero.(u[1:3])
        for (phase, x, D) in ((1, X, Dl), (2, Y, Dv))
            mean_mass = vec(sum(mw.*x; dims = 1))
            c = (rho[phase, 1]/mean_mass[1] + rho[phase, 2]/mean_mass[2])/2
            result .-= D[:, 1].*minimum(S[phase, :]).*c.*mw.*(x[:, 2] - x[:, 1])
        end
        return result
    end
    expected = reference_flux(u)
    @test flux(u) ≈ expected rtol = 1e-12
    @test flux(u, TPFA(2, 1, -1)) ≈ -expected rtol = 1e-12
    @test all(iszero, flux(u, grad, [0.0 1.0; 1.0 0.0]))
    derivatives = JutulDarcy.ForwardDiff.jacobian(flux, u)
    @test derivatives ≈ JutulDarcy.ForwardDiff.jacobian(reference_flux, u) rtol = 1e-12
    mpfa = Jutul.NFVM.NFVMLinearDiscretization(1.0; left = 1, right = 2)
    @test flux(u, mpfa) ≈ expected rtol = 1e-12
    @test JutulDarcy.ForwardDiff.jacobian(x -> flux(x, mpfa), u) ≈ derivatives rtol = 1e-12
    finite_difference = similar(derivatives)
    for component in eachindex(u)
        plus, minus = copy(u), copy(u)
        plus[component] += 1e-7
        minus[component] -= 1e-7
        finite_difference[:, component] = (reference_flux(plus) - reference_flux(minus))/(2e-7)
    end
    @test derivatives ≈ finite_difference rtol = 1e-7
end

@testset "Flash result access through AD state wrappers" begin
    x = JutulDarcy.SVector(0.3, 0.7)
    f = MultiComponentFlash.FlashedMixture2Phase(MultiComponentFlash.single_phase_l,
        x, 0.0, x, x, 1.0, 1.0)
    state = (FlashResults = [f],)
    T = typeof(Jutul.get_ad_entity_scalar(0.0, 1, 1; tag = Cells()))
    @test_throws ArgumentError Jutul.as_value(state).FlashResults[1]
    @test_throws ArgumentError local_ad(state, 1, T).FlashResults[1]
    @test_throws ArgumentError local_ad(state, 1, T, :Pressure).FlashResults[1]
end

@testset "Component diffusion parameter setup" begin
    mesh = CartesianMesh((2, 1, 1), (2.0, 1.0, 1.0))
    mixture = MultiComponentMixture(["Methane", "CarbonDioxide", "n-Decane"])
    eos = GenericCubicEOS(mixture, PengRobinson())
    system = MultiPhaseCompositionalSystemLV(eos, (LiquidPhase(), VaporPhase()))

    # A coefficient per phase and cell remains a valid input. Each phase's
    # coefficient is repeated across components, with geometric conductance
    # equal to porosity for this unit-length, unit-area face.
    domain = reservoir_domain(mesh; porosity = 0.2,
        diffusion = [1e-9 1e-9; 3e-9 3e-9])
    model, parameters = diffusion_test_model(domain, system)
    @test Jutul.values_per_entity(model, JutulDarcy.Diffusivities()) == 3
    @test parameters[:LiquidDiffusivities] ≈ fill(0.2e-9, 3, 1)
    @test parameters[:VaporDiffusivities] ≈ fill(0.6e-9, 3, 1)
    @test !haskey(model.parameters, :AqueousDiffusivities)
    @test !haskey(model.parameters, :Diffusivities)

    # Component rows and phase-specific overrides are unambiguous under a
    # prefix, including when the component and phase counts happen to match.
    domain[:liquid_diffusion, Cells()] = [1e-9 3e-9; 2e-9 2e-9; 4e-9 4e-9]
    domain[:vapor_diffusivities, Faces()] = reshape([4e-9, 5e-9, 6e-9], 3, 1)
    model, parameters = diffusion_test_model(domain, system)
    @test parameters[:LiquidDiffusivities][:, 1] ≈ [0.3e-9, 0.4e-9, 0.8e-9]
    @test parameters[:VaporDiffusivities][:, 1] == [4e-9, 5e-9, 6e-9]

    # Legacy face conductances are already geometric values, so expansion
    # must not multiply them by area/distance or porosity a second time.
    domain = reservoir_domain(mesh; porosity = 0.2)
    domain[:diffusivities, Faces()] = reshape([7e-9, 8e-9], 2, 1)
    model, parameters = diffusion_test_model(domain, system)
    @test parameters[:LiquidDiffusivities] == fill(7e-9, 3, 1)
    @test parameters[:VaporDiffusivities] == fill(8e-9, 3, 1)

    domain = reservoir_domain(mesh; porosity = 0.2, diffusion = [2e-9, 2e-9])
    model, parameters = diffusion_test_model(domain, system)
    @test parameters[:LiquidDiffusivities] ≈ fill(0.4e-9, 3, 1)
    @test parameters[:VaporDiffusivities] ≈ fill(0.4e-9, 3, 1)

    domain = reservoir_domain(mesh; porosity = 0.2,
        liquid_diffusion = 2e-9, aqueous_diffusion = 3e-9)
    model, parameters = diffusion_test_model(domain, system)
    @test parameters[:LiquidDiffusivities] ≈ fill(0.4e-9, 3, 1)
    @test !haskey(parameters, :AqueousDiffusivities)
    @test !haskey(parameters, :VaporDiffusivities)

    single_phase = SinglePhaseSystem(AqueousPhase())
    model, parameters = diffusion_test_model(domain, single_phase)
    @test parameters[:AqueousDiffusivities] ≈ fill(0.6e-9, 1, 1)

    three_phase = ImmiscibleSystem((AqueousPhase(), LiquidPhase(), VaporPhase()))
    domain = reservoir_domain(mesh; porosity = 0.2,
        aqueous_diffusion = [1e-9 1e-9; 2e-9 2e-9; 3e-9 3e-9],
        liquid_diffusion = [4e-9, 4e-9], vapor_diffusion = [5e-9, 5e-9])
    model, parameters = diffusion_test_model(domain, three_phase)
    @test parameters[:AqueousDiffusivities][:, 1] ≈ [0.2e-9, 0.4e-9, 0.6e-9]
    @test parameters[:LiquidDiffusivities] ≈ fill(0.8e-9, 3, 1)
    @test parameters[:VaporDiffusivities] ≈ fill(1e-9, 3, 1)

    # Per-phase interpretation wins for legacy matrices even when nph == ncomp.
    domain = reservoir_domain(mesh; porosity = 0.2,
        diffusion = [1e-9 1e-9; 2e-9 2e-9; 3e-9 3e-9])
    model, parameters = diffusion_test_model(domain, three_phase)
    @test parameters[:AqueousDiffusivities] ≈ fill(0.2e-9, 3, 1)
    @test parameters[:LiquidDiffusivities] ≈ fill(0.4e-9, 3, 1)
    @test parameters[:VaporDiffusivities] ≈ fill(0.6e-9, 3, 1)

    domain = reservoir_domain(mesh; porosity = 0.2,
        vapor_diffusion = [1e-9, 1e-9])
    model, parameters = diffusion_test_model(domain, system)
    @test haskey(parameters, :VaporDiffusivities)
    @test !haskey(parameters, :LiquidDiffusivities)
    @test !haskey(parameters, :AqueousDiffusivities)

    domain[:vapor_diffusion, Cells()] = fill(1e-9, 2, 2)
    @test_throws ArgumentError diffusion_test_model(domain, system)
end

@testset "Component molar diffusion flux" begin
    S = [0.2 0.4; 0.8 0.6]
    rho = [1000.0 1100.0; 5.0 7.0]
    X = [0.7 0.6; 0.2 0.2; 0.1 0.2]
    Y = [0.1 0.2; 0.6 0.5; 0.3 0.3]
    mw = [18.015, 16.043, 2.016]./1000
    Dl = reshape([1e-9, 2e-9, 3e-9], 3, 1)
    Dv = reshape([4e-9, 5e-9, 6e-9], 3, 1)
    q0 = JutulDarcy.SVector{3}(0.0, 0.0, 0.0)
    function flux(Dl, Dv, grad, saturation = S)
        liquid = mw.*X./sum(mw.*X; dims = 1)
        vapor = mw.*Y./sum(mw.*Y; dims = 1)
        q = JutulDarcy.add_phase_diffusive_component_flux(q0, Dl, 1, grad,
            saturation, rho, liquid, 1, mw, Val(3))
        return JutulDarcy.add_phase_diffusive_component_flux(q, Dv, 1, grad,
            saturation, rho, vapor, 2, mw, Val(3))
    end
    cl = (rho[1, 1]/sum(mw.*X[:, 1]) + rho[1, 2]/sum(mw.*X[:, 2]))/2
    cv = (rho[2, 1]/sum(mw.*Y[:, 1]) + rho[2, 2]/sum(mw.*Y[:, 2]))/2
    expected = -Dl[:, 1].*0.2.*cl.*mw.*(X[:, 2] - X[:, 1]) -
        Dv[:, 1].*0.6.*cv.*mw.*(Y[:, 2] - Y[:, 1])
    @test flux(Dl, Dv, TPFA(1, 2, 1)) ≈ expected
    @test flux(Dl, Dv, TPFA(2, 1, -1)) ≈ -expected
    @test all(iszero, flux(nothing, nothing, TPFA(1, 2, 1)))
    @test all(iszero, flux(nothing, Dv, TPFA(1, 2, 1), [0.2 1.0; 0.8 0.0]))
    @test all(iszero, flux(Dl, nothing, TPFA(1, 2, 1), [0.0 0.4; 1.0 0.6]))
    common = fill(1e-9, 3, 1)
    @test abs(sum(flux(common, common, TPFA(1, 2, 1))./mw)) < 1e-18

    # Black-oil Rs/Rv are surface-volume ratios: reference densities and
    # component molar masses recover the mole fractions in each phase.
    R = [0.1, 0.2]
    oil_mass, gas_mass = 0.2, 0.016
    blackoil_flux(D, grad, saturations = S) = JutulDarcy.blackoil_diffusion(
        R, saturations, rho, 800.0, 1.0, oil_mass, gas_mass,
        1, D, 1, (1, 2), grad)
    qself, qother = blackoil_flux(Dl, TPFA(1, 2, 1))
    @test qself/qother ≈ -Dl[1]*oil_mass/(Dl[2]*gas_mass)
    qself_reverse, qother_reverse = blackoil_flux(Dl, TPFA(2, 1, -1))
    @test qself_reverse ≈ -qself
    @test qother_reverse ≈ -qother
    @test all(iszero, blackoil_flux(Dl, TPFA(1, 2, 1), [0.0 0.4; 1.0 0.6]))

    qself, qother = blackoil_flux(common, TPFA(1, 2, 1))
    xleft = (800.0/oil_mass)/(800.0/oil_mass + R[1]/gas_mass)
    xright = (800.0/oil_mass)/(800.0/oil_mass + R[2]/gas_mass)
    cl = rho[1, 1]/(oil_mass*xleft + gas_mass*(1-xleft))
    cr = rho[1, 2]/(oil_mass*xright + gas_mass*(1-xright))
    expected = -1e-9*0.2*(cl + cr)/2*oil_mass*(xright - xleft)
    @test qself ≈ expected
    @test abs(qself/oil_mass + qother/gas_mass) < 1e-18
    # Do not clip derivatives at zero solution ratio.
    fraction(rs) = JutulDarcy.black_oil_phase_mole_fraction(
        800.0, 1.0, oil_mass, gas_mass, [rs], 1)
    @test JutulDarcy.ForwardDiff.derivative(fraction, 0.0) ≈ -oil_mass/(800gas_mass)
end

@testset "Black-oil diffusion molar masses" begin
    mesh = CartesianMesh((2, 1, 1), (2.0, 1.0, 1.0))
    pvt = JutulDarcy.blackoil_bench_pvt(:spe1)
    rs_max = JutulDarcy.saturated_table(pvt[:pvt][2].tab[1])
    system = StandardBlackOilSystem(rs_max = rs_max,
        phases = (AqueousPhase(), LiquidPhase(), VaporPhase()),
        reference_densities = pvt[:rhoS])
    domain = reservoir_domain(mesh)
    model, parameters = diffusion_test_model(domain, system)
    @test !haskey(parameters, :ComponentMolarMasses)
    domain[:liquid_diffusion, Cells()] = fill(1e-9, 3, 2)
    @test_throws ArgumentError diffusion_test_model(domain, system)
    domain[:component_molar_masses, nothing] = (0.018, 0.2, 0.016)
    model, parameters = diffusion_test_model(domain, system)
    @test parameters[:ComponentMolarMasses] == [0.018 0.018; 0.2 0.2; 0.016 0.016]
    @test Jutul.values_per_entity(model, JutulDarcy.ComponentMolarMasses()) == 3
    domain[:component_molar_masses, nothing] = (0.2, 0.016)
    @test_throws ArgumentError diffusion_test_model(domain, system)
    domain[:component_molar_masses, nothing] = (0.018, 0.2, 0.0)
    @test_throws ArgumentError diffusion_test_model(domain, system)
end

@testset "Black-oil molar diffusion scalar/kernel simulation" begin
    mesh = CartesianMesh((2, 1, 1), (2.0, 1.0, 1.0))
    domain = reservoir_domain(mesh; permeability = 1e-15,
        liquid_diffusion = fill(1e-5, 3, 2))
    domain[:component_molar_masses, nothing] = (0.018, 0.2, 0.016)
    pvt = JutulDarcy.blackoil_bench_pvt(:spe1)
    model, parameters = JutulDarcy.setup_reservoir_model_from_blackoil_tables(domain;
        pvtw = pvt[:pvt][1], pvto = pvt[:pvt][2], pvdg = pvt[:pvt][3],
        reference_densities = pvt[:rhoS], extra_out = true,
        extra_outputs = [:Rs, :TotalMasses, :Saturations])
    sys = reservoir_model(model).system
    pressure = 15e6
    bo = [BlackOilX(sys, pressure; sw = 0.2, rs = rs) for rs in (10.0, 20.0)]
    state0 = setup_reservoir_state(model; Pressure = pressure,
        ImmiscibleSaturation = 0.2, BlackOilUnknown = bo)
    reference = simulate_reservoir(state0, model, [86400.0]; parameters = parameters,
        info_level = -1, linear_solver = nothing, tol_cnv = 1e-8)
    kernel = simulate_reservoir(state0, model, [86400.0]; parameters = parameters,
        mode = :ka, info_level = -1, linear_solver = nothing, tol_cnv = 1e-8)
    @test length(reference.states) == length(kernel.states) == 1
    for key in (:Pressure, :Rs, :Saturations, :TotalMasses)
        @test kernel.states[end][key] ≈ reference.states[end][key] rtol = 1e-8
    end
    rs = reference.states[end][:Rs]
    @test rs[1] > 10.0
    @test rs[2] < 20.0
    @test all(isfinite, reference.states[end][:TotalMasses])
end
