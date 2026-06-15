
using GeoEnergyIO, Jutul, JutulDarcy, Test

@testset "reservoir_mesh" begin
    m1 = reservoir_mesh((10, 20, 1), (100.0, 200.0, 1.0))
    m2 = reservoir_mesh(Val(:cartesian), (10, 20, 1), (100.0, 200.0, 1.0))
    m3 = reservoir_mesh("cartesian", (10, 20, 1), (100.0, 200.0, 1.0))
    m4 = reservoir_mesh(nx = 10, ny = 20, nz = 1, Lx = 100.0, Ly = 200.0, Lz = 1.0)
    m5 = reservoir_mesh(nx = 10, ny = 20, Lx = 100.0, Ly = 200.0)
    for m in (m1, m2, m3, m4, m5)
        @test number_of_cells(m) == 200
        @test grid_dims_ijk(m) == (10, 20, 1)
        @test m isa UnstructuredMesh
        @test m.z_is_depth
        geo = tpfv_geometry(m)

        x = geo.boundary_centroids[1, :]
        y = geo.boundary_centroids[2, :]
        z = geo.boundary_centroids[3, :]

        @test maximum(z) == 1.0
        @test minimum(z) == 0.0
        @test maximum(x) == 100.0
        @test minimum(x) == 0.0
        @test maximum(y) == 200.0
        @test minimum(y) == 0.0
    end
    actnum = ones(Int, 10, 20, 1)
    actnum[1:5, 1:10, 1] .= 0
    m6 = reservoir_mesh((10, 20, 1), (100.0, 200.0, 1.0), actnum = actnum)
    m7 = reservoir_mesh((10, 20, 1), (100.0, 200.0, 1.0), actnum = vec(actnum))
    for m in (m6, m7)
        @test number_of_cells(m) == 150
        @test grid_dims_ijk(m) == (10, 20, 1)
        @test m isa UnstructuredMesh
        @test m.z_is_depth
    end

    spe1 = parse_data_file(JutulDarcy.GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA"))
    spe9 = parse_data_file(JutulDarcy.GeoEnergyIO.test_input_file_path("SPE9", "SPE9.DATA"))
    for d in (spe1, spe9)
        m = reservoir_mesh(d["GRID"])
        @test m isa UnstructuredMesh
        @test m.z_is_depth
    end

    case = setup_case_from_data_file(spe1)
    rmesh = reservoir_mesh(case)
    @test rmesh isa UnstructuredMesh
    @test reservoir_mesh(reservoir_model(case)) == rmesh
    @test reservoir_mesh(reservoir_domain(case.model)) == rmesh
    @test reservoir_mesh(case.model) == rmesh
end

