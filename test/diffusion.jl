using Jutul, JutulDarcy, MultiComponentFlash, Test

function diffusion_test_model(domain, system)
    model = setup_reservoir_model(domain, system)
    parameters = setup_parameters(model)
    return reservoir_model(model), parameters[:Reservoir]
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

@testset "Component mass diffusion flux" begin
    S = [0.2 0.4; 0.8 0.6]
    rho = [1000.0 1100.0; 5.0 7.0]
    X = [0.7 0.6; 0.2 0.2; 0.1 0.2]
    Y = [0.1 0.2; 0.6 0.5; 0.3 0.3]
    Dl = reshape([1e-9, 2e-9, 3e-9], 3, 1)
    Dv = reshape([4e-9, 5e-9, 6e-9], 3, 1)
    q0 = JutulDarcy.SVector{3}(0.0, 0.0, 0.0)
    flux(D, grad, saturation = S) = JutulDarcy.add_diffusive_component_flux(
        q0, saturation, rho, X, Y, D, 1, grad, (1, 2), Val(3))
    expected = -Dl[:, 1] .* 0.2 .* 1050.0 .* (X[:, 2] - X[:, 1]) -
        Dv[:, 1] .* 0.6 .* 6.0 .* (Y[:, 2] - Y[:, 1])
    @test flux((Dl, Dv), TPFA(1, 2, 1)) ≈ expected
    @test flux((Dl, Dv), TPFA(2, 1, -1)) ≈ -expected
    @test all(iszero, flux((nothing, nothing), TPFA(1, 2, 1)))
    @test all(iszero, flux((nothing, Dv), TPFA(1, 2, 1), [0.2 1.0; 0.8 0.0]))
    @test all(iszero, flux((Dl, nothing), TPFA(1, 2, 1), [0.0 0.4; 1.0 0.6]))
    common = fill(1e-9, 3, 1)
    @test abs(sum(flux((common, common), TPFA(1, 2, 1)))) < 1e-20

    # The black-oil path uses the two transported components of each phase.
    R = [0.1, 0.2]
    qself, qother = JutulDarcy.blackoil_diffusion(R, S, rho,
        800.0, 1.0, 1, Dl, 1, (1, 2), TPFA(1, 2, 1), SPU(1, 2))
    @test qself/qother ≈ -Dl[1]/Dl[2]
    qself_reverse, qother_reverse = JutulDarcy.blackoil_diffusion(R, S, rho,
        800.0, 1.0, 1, Dl, 1, (1, 2), TPFA(2, 1, -1), SPU(2, 1))
    @test qself_reverse ≈ -qself
    @test qother_reverse ≈ -qother

    # Expanding a legacy phase coefficient preserves the black-oil mass law.
    qself, qother = JutulDarcy.blackoil_diffusion(R, S, rho,
        800.0, 1.0, 1, common, 1, (1, 2), TPFA(1, 2, 1), SPU(1, 2))
    xleft = 800.0/(800.0 + R[1])
    xright = 800.0/(800.0 + R[2])
    expected = 1e-9*(rho[1, 1]*S[1, 1] + rho[1, 2]*S[1, 2])/2*(xleft - xright)
    @test qself ≈ expected
    @test qother ≈ -expected
end
