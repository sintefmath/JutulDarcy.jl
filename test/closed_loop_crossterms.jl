using Jutul, JutulDarcy, Test

@testset "Deprecated closed-loop cross terms" begin
    # Closed loop as separate supply and return wells of type :closed_loop,
    # coupled by cross terms (deprecated). The supply pipe flows upwards from
    # the bottom cell, and the return well is at the top cell.
    nz = 10
    domain = reservoir_domain(CartesianMesh((1, 1, nz), (1000.0, 1000.0, 100.0));
        rock_density = 2650.0, rock_heat_capacity = 900.0,
        rock_thermal_conductivity = 2.5, fluid_thermal_conductivity = 0.6,
        component_heat_capacity = 4184.0, porosity = 0.01, permeability = 1e-15)
    cells = collect(nz:-1:1)
    supply = setup_well(domain, cells;
        name = :CL_supply, type = :closed_loop, simple_well = false,
        neighborship = vcat((1:nz-1)', (2:nz)'),
        perforation_cells_well = collect(1:nz),
        well_cell_centers = domain[:cell_centroids][:, cells],
        end_nodes = [nz], WI = 0.0, radius = 20e-3)
    ret = setup_well(domain, cells[end];
        name = :CL_return, type = :closed_loop, simple_well = false,
        WI = 0.0, WIth = 0.0)
    sys = SinglePhaseSystem(AqueousPhase(), reference_density = 1000.0)
    model = @test_logs (:warn, r"deprecated") match_mode = :any setup_reservoir_model(
        domain, sys; wells = [supply, ret], thermal = true, can_shut_wells = false)
    cts = [ct.cross_term for ct in model.cross_terms]
    @test count(ct -> ct isa JutulDarcy.ClosedLoopSupplyToReturnMassCT, cts) == 1
    @test count(ct -> ct isa JutulDarcy.ClosedLoopSupplyToReturnEnergyCT, cts) == 1

    rate = 20.0si_unit(:meter)^3/si_unit(:day)
    control = Dict(
        :CL_supply => InjectorControl(TotalRateTarget(rate), [1.0],
            density = 1000.0, temperature = convert_to_si(80.0, :Celsius)),
        :CL_return => ProducerControl(BottomHolePressureTarget(5.0si_unit(:atm)))
    )
    forces = setup_reservoir_forces(model, control = control)
    state0 = setup_reservoir_state(model,
        Pressure = 5.0si_unit(:atm), Temperature = convert_to_si(10.0, :Celsius))
    res = simulate_reservoir(state0, model, fill(si_unit(:day), 2);
        forces = forces, info_level = -1)
    @test length(res.result.states) == 2
    q_in = res.wells[:CL_supply][:mass_rate][end]
    q_out = res.wells[:CL_return][:mass_rate][end]
    @test q_in > 0
    @test q_in ≈ -q_out rtol = 1e-6
end
