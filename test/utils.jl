using Jutul, JutulDarcy, Test

g = CartesianMesh((1, 1, 1), (10, 3.14, 2.71828))
r = 0.2;

@testset "Peaceman well index" begin
    @testset "Full tensor" begin
        K = [1, 0.2, 3, 5, 7, 0.15]
        @test compute_peaceman_index(g, K, r, (1, 1, 1), :x) ≈ 45.483053742963861
        @test compute_peaceman_index(g, K, r, (1, 1, 1), :y) ≈ 4.887480697644804
        @test compute_peaceman_index(g, K, r, (1, 1, 1), :z) ≈ 16.756382219371336
    end
    @testset "Diagonal" begin
        K = [0.2, 0.71, 0.314]
        @test compute_peaceman_index(g, K, r, (1, 1, 1), :x) ≈ 28.026376954943252
        @test compute_peaceman_index(g, K, r, (1, 1, 1), :y) ≈ 2.382332435163685
        @test compute_peaceman_index(g, K, r, (1, 1, 1), :z) ≈ 2.890000782901316
    end
    @testset "Scalar" begin
        K = [3.14]
        @test compute_peaceman_index(g, K, r, (1, 1, 1), :x) ≈ 1.848702389098454e+02
        @test compute_peaceman_index(g, K, r, (1, 1, 1), :y) ≈ 31.263184953683595
        @test compute_peaceman_index(g, K, r, (1, 1, 1), :z) ≈ 26.90991906669857
    end
end


import JutulDarcy: current_phase_index
@testset "current_phase_index" begin
    depths = (1.0, 2.0) # W O G
    @test current_phase_index(0.1, depths, reverse = false) == 1
    @test current_phase_index(1.01, depths, reverse = false) == 2
    @test current_phase_index(1.5, depths, reverse = false) == 2
    @test current_phase_index(2.5, depths, reverse = false) == 3

    depths = (2.0, 1.0) # WO OG
    @test current_phase_index(0.1, depths, reverse = true) == 3 # Gas
    @test current_phase_index(1.01, depths, reverse = true) == 2 # Oil
    @test current_phase_index(1.5, depths, reverse = true) == 2 # Oil
    @test current_phase_index(2.5, depths, reverse = true) == 1 # Water
end

@testset "nnc" begin
    l = [1, 2, 3]
    r = [2, 3, 4]
    N = [l'; r']
    neighbors = [(l[i], r[i]) for i in eachindex(l)]
    T = [10.0, 5.0, 3.0]
    T_t = [0.1, 0.2, 0.3]

    m = UnstructuredMesh(CartesianMesh((5, 5, 1)))
    nnc1 = setup_nnc_connections(m, l, r, T, T_t)
    nnc2 = setup_nnc_connections(m, N, T, T_t)
    @test nnc1.cells == nnc2.cells
    @test nnc1.trans_flow == nnc2.trans_flow
    @test nnc1.trans_thermal == nnc2.trans_thermal

    nnc3 = setup_nnc_connections(m, neighbors, T, T_t)
    @test nnc1.cells == nnc3.cells
    @test nnc1.trans_flow == nnc3.trans_flow
    @test nnc1.trans_thermal == nnc3.trans_thermal

    nnc4 = setup_nnc_connections(m, l, r, T)
    @test nnc1.cells == nnc4.cells
    @test all(isequal(0.0), nnc4.trans_thermal)

    @test_throws "ArgumentError: Left neighbor index 26 at connection 1 exceeds number of cells 25" JutulDarcy.setup_nnc_connections(m, [26], [1], T, T_t)
    @test_throws "ArgumentError: Right neighbor index 0 at connection 1 must be positive" JutulDarcy.setup_nnc_connections(m, [1], [0], T, T_t)

    d0 = reservoir_domain(m)
    d = reservoir_domain(m, nnc = nnc1)

    @test number_of_cells(d) == number_of_cells(d)
    @test number_of_faces(d) == number_of_faces(d0) + 3
    @test number_of_half_faces(d) == number_of_half_faces(d0) + 6

    T_computed0 = JutulDarcy.reservoir_transmissibility(d0)
    T_computed = JutulDarcy.reservoir_transmissibility(d)
    @test length(T_computed) == length(T_computed0) + 3
    @test all(T_computed[1:end-3] .== T_computed0)
    @test T_computed[end-2] == 10.0
    @test T_computed[end-1] == 5.0
    @test T_computed[end] == 3.0

    Tt_computed0 = JutulDarcy.reservoir_conductivity(d0)
    Tt_computed = JutulDarcy.reservoir_conductivity(d)
    @test length(Tt_computed) == length(Tt_computed0) + 3
    @test all(Tt_computed[1:end-3] .== Tt_computed0)
    @test Tt_computed[end-2] == 0.1
    @test Tt_computed[end-1] == 0.2
    @test Tt_computed[end] == 0.3
end
