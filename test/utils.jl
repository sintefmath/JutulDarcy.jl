using Jutul, JutulDarcy

g = CartesianMesh((1, 1, 1), (10, 3.14, 2.71828))
r = 0.2;

@testset "Peaceman well index" begin
    @testset "Full tensor" begin
        K = [1, 0.2, 3, 5, 7, 0.15]
        @test compute_peaceman_index(g, K, r, (1, 1, 1), dir = :x) ≈ 45.483053742963861
        @test compute_peaceman_index(g, K, r, (1, 1, 1), dir = :y) ≈ 4.887480697644804
        @test compute_peaceman_index(g, K, r, (1, 1, 1), dir = :z) ≈ 16.756382219371336
    end
    @testset "Diagonal" begin
        K = [0.2, 0.71, 0.314]
        @test compute_peaceman_index(g, K, r, (1, 1, 1), dir = :x) ≈ 28.026376954943252
        @test compute_peaceman_index(g, K, r, (1, 1, 1), dir = :y) ≈ 2.382332435163685
        @test compute_peaceman_index(g, K, r, (1, 1, 1), dir = :z) ≈ 2.890000782901316
    end
    @testset "Scalar" begin
        K = [3.14]
        @test compute_peaceman_index(g, K, r, (1, 1, 1), dir = :x) ≈ 1.848702389098454e+02
        @test compute_peaceman_index(g, K, r, (1, 1, 1), dir = :y) ≈ 31.263184953683595
        @test compute_peaceman_index(g, K, r, (1, 1, 1), dir = :z) ≈ 26.90991906669857
    end
end

