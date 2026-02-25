using Jutul, JutulDarcy, GeoEnergyIO, Test

@testset "Sensitivity of equilibrium region" begin

    nx = 1
    nz = 30
    Darcy, bar, yr, Kelvin = si_units(:darcy, :bar, :year, :Kelvin)
    cmesh = CartesianMesh((nx, 1, nz), (20.0, 10.0, 1000.0))
    cmesh = UnstructuredMesh(cmesh, z_is_depth = true)
    for (i, pt) in enumerate(cmesh.node_points)
        cmesh.node_points[i] = pt .+ [0.0, 0.0, 800.0]
    end
    reservoir = reservoir_domain(cmesh, porosity = 0.25, permeability = 0.1*Darcy)
    ncells = number_of_cells(cmesh)

    W = AqueousPhase()
    O = LiquidPhase()
    G = VaporPhase()
    rhoS = [1000.0, 700.0, 100.0]

    sys_3ph = ImmiscibleSystem((W, O, G), reference_densities = rhoS)
    rho = ConstantCompressibilityDensities(sys_3ph, 1.5bar, rhoS, 1e-5/bar)
    model_3ph = setup_reservoir_model(reservoir, sys_3ph)
    set_secondary_variables!(model_3ph[:Reservoir], PhaseMassDensities = rho)


    # eql = EquilibriumRegion(model_3ph, 100*bar, 1000.0, woc = 1600.0, goc = 1100.0)
    # state0_3_ph_eql = setup_reservoir_state(model_3ph, eql)

    spe9 = GeoEnergyIO.test_input_file_path("SPE9", "SPE9.DATA");
    case_spe9 = setup_case_from_data_file(spe9)
    model_bo = setup_reservoir_model(reservoir, case_spe9.model);

    function mean_pressure(m, s, dt, step_no, forces)
        return sum(s[:Reservoir][:Pressure])/length(s[:Reservoir][:Pressure])/si_unit(:bar)
    end

    function gas_saturation(m, s, dt, step_no, forces)
        return sum(s[:Reservoir][:Saturations][3, :])
    end

    function setup_case(prm, step_info = missing; variant = "3ph")
        if variant == "3ph"
            model = model_3ph
        else
            model = model_bo
        end

        dt = [1000.0]
        p_datum = prm[:p_datum]
        datum_depth = prm[:datum_depth]
        kwarg = Dict()
        for (k, v) in pairs(prm)
            if k == :p_datum || k == :datum_depth
                continue
            else
                kwarg[k] = v
            end
        end
        eql = EquilibriumRegion(model, p_datum, datum_depth; kwarg...)
        s = setup_reservoir_state(model, eql)
        return JutulCase(model, dt, state0 = s)
    end

    prm0 = Dict(
        :p_datum => 100*bar,
        :datum_depth => 1000.0,
        :woc => 1600.0,
        :goc => 1100.0,
        :T_z => (z) -> 273.15 .+ 0.01 .* z, # Will still test stuff even if models are not thermal
        :rs => 15.0
    )

    f_bo = (arg...) -> setup_case(arg..., variant = "bo")
    f_3ph = (arg...) -> setup_case(arg..., variant = "3ph")

    function test_gradient(F, obj)
        dopt = Jutul.DictOptimization.DictParameters(prm0, F, strict = false)
        free_optimization_parameters!(dopt)
        problem = Jutul.DictOptimization.JutulOptimizationProblem(dopt, obj, info_level = -1);

        o, dodx = problem()
        dodx_num = similar(dodx)
        for i in eachindex(dodx)
            dodx_num[i] = Jutul.DictOptimization.finite_difference_gradient_entry(problem, index = i)
        end
        for i in eachindex(dodx)
            @test dodx[i] ≈ dodx_num[i] atol = 1e-6
        end
    end

    test_gradient(f_3ph, mean_pressure)
    test_gradient(f_bo, mean_pressure)

    test_gradient(f_3ph, gas_saturation)
    test_gradient(f_bo, gas_saturation)
end
