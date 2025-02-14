using Jutul, JutulDarcy, Test, Statistics

function solve_btes(;
        nc = 3,
        btes_type = :simple,
        use_tables = true,
        temp_charge = convert_to_si(50.0, :Celsius),
        temp_discharge = convert_to_si(10.0, :Celsius),
        n_step = 10
    )

    darcy, atm, litre, year, second = si_units(:darcy, :atm, :litre, :year, :second)

    g = CartesianMesh((nc, nc, nc), (50.0, 50.0, 50.0))
    d = reservoir_domain(g,
        permeability = 10e-3darcy,
        porosity = 0.05,
        rock_thermal_conductivity = 2.0,
        fluid_thermal_conductivity = 0.6,
    )

    btes_sup, btes_ret = setup_vertical_btes_well(d, 2, 2, name = :BTES, btes_type = btes_type)
    if use_tables
        sys = :geothermal
    else
        sys = SinglePhaseSystem(AqueousPhase(), reference_density = 1000.0)
    end
    model, = setup_reservoir_model(d, sys;
        wells = [btes_sup, btes_ret], thermal = true)
    rmodel = reservoir_model(model)
    sys = rmodel.system
    rho_w = first(JutulDarcy.reference_densities(sys))

    k = map(c -> cell_ijk(g, c)[3], 1:number_of_cells(g))
    top = findall(k .== 1)
    bc = FlowBoundaryCondition.(top, 1atm, temp_discharge)

    rate_target = TotalRateTarget(1.0litre/second)
    ctrl_charge  = InjectorControl(rate_target, [1.0], density = 1000.0, temperature = temp_charge)
    ctrl_discharge = InjectorControl(rate_target, [1.0], density = 1000.0, temperature = temp_discharge)

    bhp_target = BottomHolePressureTarget(1atm)
    ctrl_prod = ProducerControl(bhp_target)

    control_charge    = Dict(:BTES_supply => ctrl_charge, :BTES_return => ctrl_prod);
    control_discharge = Dict(:BTES_supply => ctrl_discharge, :BTES_return => ctrl_prod);

    forces_charge = setup_reservoir_forces(model, control = control_charge, bc = bc)
    forces_discharge = setup_reservoir_forces(model, control = control_discharge, bc = bc)

    time = 1year
    dt = fill(time, n_step)
    forces = vcat(fill(forces_charge, n_step), fill(forces_discharge, n_step))
    dt = repeat(dt, 2)

    state0 = setup_reservoir_state(model, Pressure = 10.0*atm, Temperature = temp_discharge)

    states, reports = simulate(state0, model, dt,
        forces = forces,
        info_level = -1,
        error_on_incomplete = true,
        failure_cuts_timestep = false
    )

    return states, reports

end

@testset "geothermal btes" begin
    for use_tables in [true, false]
        for btes_type in [:simple, :u1]
            n_step = 10
            temp_charge = convert_to_si(50.0, :Celsius)
            temp_discharge = convert_to_si(10.0, :Celsius)
            F() = solve_btes(
                btes_type = btes_type, 
                temp_charge = temp_charge,
                temp_discharge = temp_discharge,
                n_step = n_step,
                use_tables = use_tables
            )

            if use_tables && btes_type == :simple
                # This combo is not yet converging.
                @test_broken F()
            else
                states, = F()
                T_res = states[n_step][:Reservoir][:Temperature]
                @test temp_discharge < maximum(T_res) <= temp_charge

                T_btes = states[n_step][:BTES_return][:Temperature]
                @test isapprox(maximum(T_btes), temp_charge, atol = 2.0)

                T_res = states[end][:Reservoir][:Temperature]
                @test temp_discharge < maximum(T_res) <= temp_charge

                T_btes = states[end][:BTES_return][:Temperature]
                @test isapprox(minimum(T_btes), temp_discharge, atol = 2.0)
            end
        end
    end
end
