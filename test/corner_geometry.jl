using Test, Jutul, GeoEnergyIO, JutulDarcy
import JutulDarcy.CornerPointGeometry as CPGD

@testset "Corner-point reservoir geometry" begin
    for rotation in (0.0,40.0), throw in (0.0,5.0)
        nx, ny, nz = 2, 1, 2
        coord = Float64[]
        for j in 0:ny, i in 0:nx
            x, y = 100.0i, 100.0j
            xr, yr = cosd(rotation)*x-sind(rotation)*y, sind(rotation)*x+cosd(rotation)*y
            append!(coord, (xr,yr,0.0,xr,yr,3000.0))
        end
        zcorn = [1000.0+(k-1+dk)*10+(i-1)*throw
            for di in 0:1, i in 1:nx, dj in 0:1, j in 1:ny, dk in 0:1, k in 1:nz]
        grid = Dict{String,Any}("cartDims"=>(nx,ny,nz), "COORD"=>coord, "ZCORN"=>vec(zcorn), "ACTNUM"=>trues(nx,ny,nz))
        dom = reservoir_domain(mesh_from_grid_section(grid); permeability=1e-13,porosity=0.2)
        CPGD.set_cpgrid_geometry!(dom,grid)
        @test all(isapprox.(dom[:volumes], 100000))
        @test all(isapprox.(JutulDarcy.peaceman_cell_dims(dom,1), (100,100,10)))
        @test JutulDarcy.compute_peaceman_index(dom,1e-13,0.1,1) ≈
            JutulDarcy.compute_peaceman_index((100.0,100.0,10.0),1e-13,0.1)
        well = setup_well(dom,[1]; simple_well=true)
        @test all(isapprox.(well[:cell_dims][1], (100,100,10)))
        for cell in (1,(1,1,1))
            well = setup_well(dom,cell;simple_well=true)
            @test all(isapprox.(well[:cell_dims][1], (100,100,10)))
        end
        # Isotropic permeability makes the fault result invariant to rotation.
        N = dom[:neighbors]
        f = only(findall(f -> Set(N[:,f]) == Set((1,2)), axes(N,2)))
        @test JutulDarcy.reservoir_transmissibility(dom)[f] ≈ 1e-13*dom[:areas][f]/100
    end
    top = [1000.0+0.2x+0.1y for x in (0.0,100.0,200.0), y in (0.0,100.0,200.0)]
    grid = cpgrid_from_horizons([0.0,100.0,200.0], [0.0,100.0,200.0], [top,top .+ 10])
    dom = reservoir_domain(mesh_from_grid_section(grid);permeability=1e-13,porosity=0.2)
    CPGD.set_cpgrid_geometry!(dom,grid)
    p = dom[:cpgrid_corners][1]
    @test p[2][3]-p[1][3] ≈ 20
    @test p[3][3]-p[1][3] ≈ 10
    for f in axes(dom[:neighbors],2)
        l,r = dom[:neighbors][:,f]
        @test dom[:face_centroids][:,f] ≈ (dom[:cell_centroids][:,l]+dom[:cell_centroids][:,r])/2
    end
    coarse = JutulDarcy.coarsen_reservoir(dom, [1,1,2,2])
    @test !haskey(coarse, :cpgrid_corners)
    @test !haskey(coarse, :half_face_centroids)
    @test coarse[:volumes] ≈ [200000,200000]
    @test only(JutulDarcy.reservoir_transmissibility(coarse)) > 0.0
end
