using Jutul, JutulDarcy
using Test

function test_multisegment_well_orientation_invariance()
    # Grid: 3x3x25 cells covering 600x600x1000 domain
    nx, ny, nz = 3, 3, 25
    grid_dims = (600.0, 600.0, 1000.0)

    day = 3600 * 24
    bar = 1e5

    g = CartesianMesh((nx, ny, nz), grid_dims)
    Darcy = 9.869232667160130e-13
    domain = reservoir_domain(g, permeability = 1.0*Darcy, porosity = 0.2)
    cell_centers = domain[:cell_centroids]

    # Producer: center column; Injector: corner column
    prod_perf_cells = [cell_index(g, (nx ÷ 2 + 1, ny ÷ 2 + 1, k)) for k in 1:nz]
    inj_perf_cells  = [cell_index(g, (1, 1, k)) for k in 1:nz]

    # Two-phase immiscible system with strong density contrast for pronounced buoyancy
    phases = (LiquidPhase(), VaporPhase())
    rhoLS = 1000.0  # Liquid (dense)
    rhoGS = 50.0    # Gas (light) - 20x lighter than liquid
    rhoS = [rhoLS, rhoGS]
    sys = ImmiscibleSystem(phases, reference_densities = rhoS)
    c = [1e-6/bar, 1e-4/bar]
    rho = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
    
    # Explicit well topology: segments 1->2->...->nz (top to bottom)
    num_perfs = nz
    seg_from = collect(1:num_perfs - 1)
    seg_to   = collect(2:num_perfs)
    connectivity_normal = vcat(seg_from', seg_to')
    perf_to_wellcell = collect(1:num_perfs)
    prod_well_centers = cell_centers[:, prod_perf_cells]
    inj_well_centers  = cell_centers[:, inj_perf_cells]

    n_closed = nz ÷ 4  # upper 25% of perforations: WI = 0

    function make_well(perf_cells, well_centers, connectivity, well_name)
        well = setup_well(
            domain, perf_cells,
            neighborship           = connectivity,
            perforation_cells_well = perf_to_wellcell,
            well_cell_centers      = well_centers,
            name        = well_name,
            simple_well = false,
            radius      = 0.05,
            roughness   = 1e-4,
        )
        well[:well_index, Perforations()][1:n_closed] .= 0.0
        return well
    end

    function make_model(prod_connectivity, inj_connectivity)
        prod_well = make_well(prod_perf_cells, prod_well_centers, prod_connectivity, :Producer)
        inj_well  = make_well(inj_perf_cells,  inj_well_centers,  inj_connectivity,  :Injector)
        # Each call independently computes SegmentConnectionGravityDifference
        # from each well's own neighborship.
        model = setup_reservoir_model(
            domain, sys,
            wells = [prod_well, inj_well]
        )
        replace_variables!(model, PhaseMassDensities = rho)
        return model
    end

    controls = Dict(
        :Producer => ProducerControl(BottomHolePressureTarget(150*bar)),
        :Injector => InjectorControl(BottomHolePressureTarget(250*bar), [1.0, 0.0], density = rhoLS),
    )
    dt = [10.0] * day

    function run_simulation(prod_connectivity, inj_connectivity)
        model = make_model(prod_connectivity, inj_connectivity)
        state0   = setup_reservoir_state(model, Pressure = 200*bar, Saturations = [0.8, 0.2])
        forces   = setup_reservoir_forces(model, control = controls)
        sim, cfg = setup_reservoir_simulator(model, state0, info_level = -1)
        states, reports = simulate!(sim, dt, forces = forces, config = cfg)
        return states, reports
    end

    # Apply orientation flips to a copy of a connectivity matrix
    function apply_flips(base, indices)
        conn = copy(base)
        for i in indices
            conn[1, i], conn[2, i] = conn[2, i], conn[1, i]
        end
        return conn
    end

    # Build a sign vector: -1 for each flipped segment, +1 otherwise
    function orientation_signs(flipped_indices)
        s = ones(num_perfs - 1)
        for i in flipped_indices
            s[i] = -1.0
        end
        return s
    end

    # Reference: both wells with normal topology
    states_ref, reports_ref = run_simulation(connectivity_normal, connectivity_normal)

    # Compare a test run against the reference.
    # Reservoir state and well node pressures must be identical.
    # TotalMassFlux[i] acquires a sign flip for each reversed segment.
    function check_invariance(states_test, prod_signs, inj_signs, label)
        @testset "$label" begin
            @test isapprox(
                states_ref[end][:Reservoir][:Pressure],
                states_test[end][:Reservoir][:Pressure],
                rtol = 1e-6, atol = 1e-2*bar,
            )
            @test isapprox(
                states_ref[end][:Reservoir][:Saturations],
                states_test[end][:Reservoir][:Saturations],
                rtol = 1e-5, atol = 1e-6,
            )
            @test isapprox(
                states_ref[end][:Producer][:Pressure],
                states_test[end][:Producer][:Pressure],
                rtol = 1e-5, atol = 1e-2*bar,
            )
            @test isapprox(
                states_ref[end][:Injector][:Pressure],
                states_test[end][:Injector][:Pressure],
                rtol = 1e-5, atol = 1e-2*bar,
            )
            # TotalMassFlux[i] > 0 means flow from N[1,i] toward N[2,i].
            # Reversing segment i gives V_test[i] = -V_ref[i].
            flux_P_ref  = states_ref[end][:Producer][:TotalMassFlux]
            flux_P_test = states_test[end][:Producer][:TotalMassFlux]
            @test isapprox(flux_P_ref, prod_signs .* flux_P_test, rtol = 1e-5, atol = 1e-10)
            flux_I_ref  = states_ref[end][:Injector][:TotalMassFlux]
            flux_I_test = states_test[end][:Injector][:TotalMassFlux]
            @test isapprox(flux_I_ref, inj_signs .* flux_I_test, rtol = 1e-5, atol = 1e-10)
        end
    end

    @testset "Multisegment well orientation invariance" begin
        # Case 1: middle portion reversed (segments nz/4 to 3*nz/4)
        mid_range = max(1, nz ÷ 4):min(nz - 1, 3 * nz ÷ 4)
        states_mid, reports_mid = run_simulation(
            apply_flips(connectivity_normal, mid_range),
            apply_flips(connectivity_normal, mid_range),
        )
        check_invariance(states_mid,
            orientation_signs(mid_range), orientation_signs(mid_range),
            "Middle portion reversed",
        )

        # Case 2: all segments reversed
        all_range = 1:num_perfs - 1
        states_all, reports_all = run_simulation(
            apply_flips(connectivity_normal, all_range),
            apply_flips(connectivity_normal, all_range),
        )
        check_invariance(states_all,
            orientation_signs(all_range), orientation_signs(all_range),
            "All segments reversed",
        )

        # Case 3: every second segment reversed
        alt_range = 2:2:num_perfs - 1
        states_alt, reports_alt = run_simulation(
            apply_flips(connectivity_normal, alt_range),
            apply_flips(connectivity_normal, alt_range),
        )
        check_invariance(states_alt,
            orientation_signs(alt_range), orientation_signs(alt_range),
            "Every second segment reversed",
        )
    end
end

@testset "Multisegment wells" begin
    test_multisegment_well_orientation_invariance()
end
