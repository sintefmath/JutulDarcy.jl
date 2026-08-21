using JutulDarcy, Jutul, Test
@testset "Well constructor" begin
    Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day);
    nx = ny = 5
    nz = 4
    dims = (nx, ny, nz)
    g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
    domain = reservoir_domain(g, permeability = 0.1*Darcy, porosity = 0.2)

    is_simple = false
    Inj = setup_well(domain, [(nx, ny, 1)], name = :Injector, reference_depth = -1.0)
    Inj_w = physical_representation(Inj)
    @test Inj isa DataDomain
    @test Inj_w isa WellDomain
    @test Inj_w.name == :Injector

    Prod = setup_vertical_well(domain, 1, 1, name = :Producer, simple_well = is_simple)
    Prod_shifted = setup_vertical_well(domain, 1, 1, name = :Producer, simple_well = is_simple, reference_depth = -1.0)

    phases = (LiquidPhase(), VaporPhase())
    rhoLS = 1000.0
    rhoGS = 100.0
    rhoS = [rhoLS, rhoGS] .* kg/meter^3
    sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    model, prm = setup_reservoir_model(domain, sys, wells = [Inj, Prod], thermal = true, extra_out = true);

    @test prm[:Injector][:FluidVolume][1] ≈ 0.392699 atol = 0.001
    @test prm[:Injector][:WellIndicesThermal][1] ≈ 27.3211 atol = 0.001
    @test prm[:Injector][:PerforationGravityDifference][1] ≈ 71.0982 atol = 0.001
    @test prm[:Injector][:WellIndices][1] ≈ 1.18321e-12 rtol = 1e-3
end

@testset "Well convenience forces" begin
    i1 = setup_injector_control(1.0, :rate, [0.0, 1.0], density = 100.0)
    @test i1 isa InjectorControl
    @test i1.target isa TotalRateTarget
    @test i1.target.value == 1.0
    @test i1.mixture_density == 100.0
    @test i1.injection_mixture == [0.0, 1.0]

    i1 = setup_injector_control(1.0, "rate", [0.0, 1.0], density = 100.0)
    @test i1 isa InjectorControl
    @test i1.target isa TotalRateTarget
    @test i1.target.value == 1.0
    @test i1.mixture_density == 100.0

    p = setup_producer_control(100.0, :bhp)
    @test p isa ProducerControl
    @test p.target isa BottomHolePressureTarget
    @test p.target.value == 100.0

    p = setup_producer_control(100.0, "wbhp")
    @test p isa ProducerControl
    @test p.target isa BottomHolePressureTarget
    @test p.target.value == 100.0
end
