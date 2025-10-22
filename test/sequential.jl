import JutulDarcy.Sequential:
    phase_potential_r_index,
    phase_potential_upwind_fixed_flux,
    phase_potential_upwind_potential_differences,
    sort_tuple_indices
using Test, Jutul, JutulDarcy, HYPRE, LinearAlgebra

@testset "sequential_upwind" begin
    @testset "sort_tuple" begin
        @test sort_tuple_indices((10.0, 5.0)) == (2, 1)
        @test sort_tuple_indices((5.0, 10.0)) == (1, 2)

        tups = [
            (1, 2, 3),
            (1, 3, 2),
            (2, 1, 3),
            (2, 3, 1),
            (3, 1, 2),
            (3, 2, 1)
        ]

        vals = [10.0, 50.0, 500.0]
        for ix in tups
            V = (vals[ix[1]], vals[ix[2]], vals[ix[3]])
            ref = sort(eachindex(V), by = i -> V[i])
            @test sort_tuple_indices(V) == Tuple(ref)
        end
    end

    @testset "phase_potential_r_index" begin
        @test phase_potential_r_index(1, 1) == 0
        @test phase_potential_r_index(0, 1) == 1
        @test phase_potential_r_index(-1, 1) == 1
        @test phase_potential_r_index(-1, -1) == 2
        @test phase_potential_r_index(-1, -1) == 2
    end

    simple_upwind(l, r, flag) = ifelse(flag, r, l)

    function test_flux(q, K, g_ph, mob_l, mob_r)
        flags = phase_potential_upwind_fixed_flux(q, K, g_ph, mob_l, mob_r)
        N = length(g_ph)
        vals = Float64[]
        mobT = 0.0
        upw(flag, ph) = ifelse(flag, mob_r[ph], mob_l[ph])
        for l in 1:N
            val = q
            mobT += upw(flags[l], l)
            for j in 1:N
                if l != j
                    mob = upw(flags[j], j)
                    val += K*(g_ph[l] - g_ph[j])*mob
                end
            end
            push!(vals, val)
        end
        vals2 = phase_potential_upwind_potential_differences(q, K, g_ph, mob_l, mob_r)

        # vals = vals2
        # @test vals == vals2
        for (i, val) in enumerate(vals)
            if val > 0
                @test flags[i] == false
            elseif val < 0
                @test flags[i] == true
            end
        end
    end

    # Also make sure that these two are precompiled before the next test.
    A = (0.0, 111.11111111111111)
    B = (111.11111111111111, 0.0)
    G = (-1.0, -10.0)
    test_flux(-1.0, 1, G, A, B)
    # @test 0 == @allocated phase_potential_upwind_fixed_flux(-1.0, 1, G, A, B)
    A = (0.0, 111.0, 1.0)
    B = (10.0, 0.0, 1.0)
    G = (-1.0, -10.0, -20.0)
    test_flux(-1.0, 1, G, A, B)
    # @test 0 == @allocated phase_potential_upwind_fixed_flux(-1.0, 1, G, A, B)

    @testset "2ph" begin
        Nmob = 10
        for sgn in [-1, 1]
            for fsgn in [-1, 1]
                for lmob1 in range(0, 1000, length = Nmob)
                    for rmob1 in range(0, 1000, length = Nmob)
                        for lmob2 in range(0, 1000, length = Nmob)
                            for rmob2 in range(0, 1000, length = Nmob)
                                # 2ph
                                A = (lmob1, lmob2)
                                B = (rmob1, rmob2)
                                if !all(isequal(0), A) && !all(isequal(0), B)
                                    test_flux(0, 1, sgn.*(1.0, 2.0), A, B)
                                    test_flux(fsgn*1, 1, sgn.*(0.0, 0.0), A, B)
                                    test_flux(fsgn*0.5, 1, sgn.*(1.0, 2.0), A, B)
                                    test_flux(fsgn*1, 1, sgn.*(10.0, 1.0), A, B)
                                    test_flux(fsgn*1, 1, sgn.*(1.0, 10.0), A, B)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    ##
    @testset "3ph" begin
        Nmob = 5
        for sgn in [-1, 1]
            for fsgn in [-1, 1]
                for lmob1 in range(0, 1000, length = Nmob)
                    for rmob1 in range(0, 1000, length = Nmob)
                        for lmob2 in range(0, 1000, length = Nmob)
                            for rmob2 in range(0, 1000, length = Nmob)
                                for lmob3 in range(0, 1000, length = Nmob)
                                    for rmob3 in range(0, 1000, length = Nmob)
                                        A = (lmob1, lmob2, lmob3)
                                        B = (rmob1, rmob2, rmob3)
                                        if !all(isequal(0), A) && !all(isequal(0), B)
                                            test_flux(0, 1, sgn.*(1.0, 2.0, 3.0), A, B)
                                            test_flux(0, 1, sgn.*(1.0, 2.0, 3.0), A, B)
                                            test_flux(fsgn*100, 1, sgn.*(1.0, 2.0, 3.0), A, B)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

##
@testset "Sequential interface" begin+
    file_path = JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA")
    case = setup_case_from_data_file(file_path, 
        extra_outputs = [:PhaseMobilities, :PhaseMassDensities, :SurfaceVolumeMobilities, :Rs]
    )
    r_seq = simulate_reservoir(case,
        method = :sequential,
        target_ds = 0.1,
        info_level = -1
    );
    ix = length(r_seq.states)
    state = r_seq.states[ix]
    state = deepcopy(state)
    state = merge!(state, case.parameters[:Reservoir])
    state[:TotalSaturation] .= 1.0
    state = JutulStorage(state)
    rmodel = reservoir_model(case.model)
    nf = number_of_faces(rmodel.domain)
    neighbors = rmodel.data_domain[:neighbors]
    nph = number_of_phases(rmodel.system)

    v_total = JutulDarcy.Sequential.store_total_fluxes(rmodel, state)
    state[:TotalVolumetricFlux] = v_total
    faces = 1:nf
    flux = zeros(nph, length(faces))
    flux_s = zeros(nph, length(faces))
    for (i, face) in enumerate(faces)
        l, r = neighbors[:, face]
        q = Jutul.StaticArrays.@MVector zeros(3)
        flux_type = Jutul.DefaultFlux()
        flux_type_s = JutulDarcy.Sequential.TotalSaturationFlux()
        upw = Jutul.SPU(l, r)
        kgrad = Jutul.TPFA(l, r, 1)
        flux[:, i] = JutulDarcy.component_mass_fluxes!(q, face, state, rmodel, flux_type, kgrad, upw)
        flux_s[:, i] = JutulDarcy.component_mass_fluxes!(q, face, state, rmodel, flux_type_s, kgrad, upw)
    end

    for ph in 1:nph
        f = flux[ph, :]
        err_scaled = abs.(f - flux_s[ph, :])/norm(f, 2)
        err = norm(err_scaled, 2)
        @test err < 1e-14
        # ix = findmax(err_scaled)[2]
        # println("Phase $ph: rel err = $err (worst value: $(f[ix]) vs $(flux_s[ph, ix]) at face $ix)")
    end
end
