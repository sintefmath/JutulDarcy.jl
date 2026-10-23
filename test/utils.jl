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

@testset "Plot unit conversion factors" begin
    @test JutulDarcy.unit_system_unit_conversion("Metric", "", :volume) == (1.0, "m³")
    @test first(JutulDarcy.unit_system_unit_conversion("Field", "", :volume)) ≈ 1 / 0.3048^3
    @test first(JutulDarcy.unit_system_unit_conversion("Field", "", :gas_volume_surface)) ≈ 1 / (1000 * 0.3048^3)
end

@testset "Reservoir plot axis proportions" begin
    # Axis3 must display the requested z-to-x ratio regardless of reservoir
    # thickness. A coordinate scale (e.g. 10 for a shallow grid) is not that ratio.
    for thickness in (0.01, 50.0, 1500.0, 15000.0)
        ranges = (1500.0, 800.0, thickness)
        fixed = JutulDarcy.reservoir_plot_aspect(ranges, 1/3; use_scaling = false)
        @test fixed[3]/fixed[1] ≈ 1/3
        @test fixed[2]/fixed[1] ≈ 800/1500
    end
    @test JutulDarcy.reservoir_plot_aspect((20.0, 40.0, 2.0), 0.5;
        use_scaling = false) == (1.0, 2.0, 0.5)

    # Only scale depth in the fancy viewer; keep its existing 2D z scale.
    @test JutulDarcy.reservoir_plot_aspect((1500.0, 800.0, 50.0), 1/3;
        use_scaling = true) == (1.0, 1.0, 10.0)
    @test JutulDarcy.reservoir_plot_aspect((20.0, 40.0), 1/3;
        use_scaling = true) == (1.0, 1.0, 1.0)
    @test JutulDarcy.reservoir_plot_aspect((20.0, 40.0), 1/3;
        use_scaling = false) == (1.0, 2.0, 1.0)
end

@testset "Fancy reservoir horizontal proportions" begin
    for (x, y) in ((1000.0, 500.0), (500.0, 1000.0), (1000.0, 1000.0),
                  (1000.0, 0.1), (0.1, 1000.0))
        for ranges in ((x, y), (x, y, 100.0))
            scale = JutulDarcy.reservoir_plot_aspect(ranges, 1/3; use_scaling = true)
            # Coordinate scaling must preserve the mesh's horizontal shape.
            @test (y*scale[2])/(x*scale[1]) ≈ y/x
        end
    end
end

@testset "partitioner_input with faults" begin
    g = UnstructuredMesh(CartesianMesh((4, 4, 1)))
    Jutul.set_mesh_entity_tag!(g, Faces(), :faults, :f1, [1, 2])
    d = reservoir_domain(g)
    nf = number_of_faces(g)
    for conn in (:unit, :trans, :logtrans)
        N, T, = JutulDarcy.partitioner_input(d, conn = conn)
        @test size(N, 2) == length(T) == nf - 2
        N, T, = JutulDarcy.partitioner_input(d, conn = conn, preserve_faults = false)
        @test size(N, 2) == length(T) == nf
    end
end
