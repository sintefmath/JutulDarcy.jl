using Jutul, JutulDarcy, GeoEnergyIO, Test

@testset "State0 sensitivities" begin
    data_dir = GeoEnergyIO.test_input_file_path("SPE1")
    data_pth = joinpath(data_dir, "SPE1.DATA")
    case = setup_case_from_data_file(data_pth)
    case = case[1:10]
    case = coarsen_reservoir_case(case, (4, 3, 1), method = :ijk);

    function obj_func(m, s, dt, step_no, forces)
        return sum(s[:Reservoir][:Pressure])/length(s[:Reservoir][:Pressure])/si_unit(:bar)
    end

    cfg = optimization_config(case,
        use_scaling = true,
        rel_min = 0.01,
        rel_max = 100,
        include_state0 = true
    )

    for (mkey, mcfg) in pairs(cfg)
        for (ki, vi) in pairs(mcfg)
            if ki == :Pressure && mkey == :Reservoir
                vi[:active] = true
            else
                vi[:active] = false
            end
        end
    end

    _, sim_cfg = setup_reservoir_simulator(case;
        rtol = 1e-5,
        tol_cnv = 1e-5,
        output_substates = true
    );

    opt_setup = setup_parameter_optimization(case, obj_func, cfg;
        print = 0,
        param_obj = true,
        config = sim_cfg
    );

    x = deepcopy(opt_setup.x0)
    F0 = opt_setup.F!(x)

    # Adjoint derivatives
    dF_adj = opt_setup.dF!(similar(x), x)

    # Numerical derivatives
    ϵ = 1e-6
    dF_num = similar(x)
    for i in eachindex(x)
        x_d = copy(x)

        x_d[i] += ϵ
        F_d = opt_setup.F!(x_d)

        x_d[i] -= 2*ϵ
        F_d2 = opt_setup.F!(x_d)
        dF_num[i] = (F_d - F_d2)/(2*ϵ)
    end

    for i in eachindex(dF_adj, dF_num)
        @test isapprox(dF_num[i], dF_adj[i], rtol = 0.05)
    end

end