
import JutulDarcy.Sequential:
    phase_potential_r_index,
    phase_potential_upwind_fixed_flux,
    phase_potential_upwind_potential_differences,
    sort_tuple
using Test
@testset "sequential_upwind" begin
    @testset "sort_tuple" begin
        @test sort_tuple((10.0, 5.0)) == (2, 1)
        @test sort_tuple((5.0, 10.0)) == (1, 2)

        combos = [
            ((1, 2, 3), (1, 2, 3)),
            ((1, 3, 2), (1, 3, 2)),
            ((2, 1, 3), (2, 1, 3)),
            ((2, 3, 1), (3, 1, 2)),
            ((3, 1, 2), (2, 3, 1)),
            ((3, 2, 1), (3, 2, 1))
        ]

        vals = [10.0, 50.0, 500.0]
        for (ix, ref) in combos
            @test sort_tuple((vals[ix[1]], vals[ix[2]], vals[ix[3]])) == ref
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

    function test_flux(q, K, g_ph, mob_l, mob_r; print = false)
        flags = phase_potential_upwind_fixed_flux(q, K, g_ph, mob_l, mob_r, print)
        N = length(g_ph)
        vals = Float64[]
        for i in 1:N
            val = q
            for j in 1:N
                if i != j
                    if flags[j]
                        mob = mob_l[j]
                    else
                        mob = mob_r[j]
                    end
                    val += K*(g_ph[i] - g_ph[j])*mob
                end
            end
            push!(vals, val)
        end
        vals2 = phase_potential_upwind_potential_differences(q, K, g_ph, mob_l, mob_r)

        vals = vals2
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
