using Test, JutulDarcy, Jutul

function test_kr_bounds(relperm, n)
    for ph in 1:n
        S = zeros(n, 1)
        S[ph] = 1.0

        kr = similar(S)
        JutulDarcy.update_kr!(kr, relperm, nothing, S, entity_eachindex(kr))
        for i in 1:n
            if i == ph
                @test kr[i] > 0
            else
                @test kr[i] == 0
            end
        end
    end
end

function kr_test_sat(n)
    if n == 2
        S = zeros(2, 3)
        S[1, 1] = 1.0
        S[1, 2] = 0.5
        S[1, 3] = 0.1
        for i = 1:3
            S[2, i] = 1 - S[1, i]
        end
    else
        error("Bad input $n")
    end
    return S
end

function test_brooks_corey_kr()
    bc = BrooksCoreyRelativePermeabilities(2, [2.0, 3.0], [0.2, 0.3], [0.9, 1.0])
    S = kr_test_sat(2)
    kr = similar(S)
    @test JutulDarcy.update_kr!(kr, bc, nothing, S, entity_eachindex(kr)) ≈ [0.9 0.324 0.0; 0.0 0.064 1.0]
    test_kr_bounds(bc, 2)
end

function test_standard_kr()
    S = kr_test_sat(2)
    kr = similar(S)

    kr_1 = S -> S^2
    kr_2 = S -> S
    rel_single = JutulDarcy.RelativePermeabilities((kr_1, kr_2))
    JutulDarcy.update_kr!(kr, rel_single, nothing, S, entity_eachindex(kr))
    for i in axes(kr, 2)
        @test kr[1, i] == kr_1(S[1, i])
        @test kr[2, i] == kr_2(S[2, i])
    end

    # Quadratic in first region, linear in second
    rel_regs = JutulDarcy.RelativePermeabilities(((kr_1, kr_2), (kr_1, kr_2)), regions = [1, 2])
    S = repeat([0.5], 2, 2)
    kr = similar(S)
    JutulDarcy.update_kr!(kr, rel_regs, nothing, S, entity_eachindex(kr))
    @test kr[1, 1] == kr[2, 1] ≈ 0.25
    @test kr[1, 2] == kr[2, 2] ≈ 0.5
end

function test_rel_perm_wrapper()
    s = [0.1, 0.15, 0.2, 0.8, 1.0]
    kr = [0.0, 0.0, 0.4, 0.9, 0.9]
    relperm = PhaseRelativePermeability(s, kr)
    @testset "Detection of points" begin
        @test relperm.connate == 0.1
        @test relperm.k_max == 0.9
        @test relperm.s_max == 0.8
        @test relperm.critical == 0.15
    end
    kro = @. (1 - s)^2
    swof = hcat(s, kr, kro)
    krw, krow = table_to_relperm(swof)
    @testset "SWOF-like table" begin
        @test krw.(s) ≈ swof[:, 2]
        @test krow.(1 .- s) ≈ swof[:, 3]
    end
end

@testset "RelativePermeabilities" begin
    @testset "BrooksCorey" begin
        test_brooks_corey_kr()
    end
    @testset "Standard kr" begin
        test_standard_kr()
    end
    @testset "Wrapper" begin
        test_rel_perm_wrapper()
    end
end
