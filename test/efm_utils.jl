@testset "Shared EFM and DFM property setup" begin
    matrix = reservoir_domain(reservoir_mesh((2, 2, 2), (2.0, 2.0, 2.0));
        permeability = repeat([2.0, 3.0, 4.0, 0.1, 0.2, 0.3].*1e-12, 1, 8))
    matrix[:region, Cells()] = collect(1:8)
    matrix[:metadata, nothing] = Dict(:values => [1, 2])
    objects = AbstractEmbeddedFracture[
        RectangleFracture([1, 1, 1], [1, 0, 0], (2, 2); aperture = 1e-3),
        RectangleFracture([1, 1, 1], [0, 1, 0], (2, 2); aperture = 2e-3)]
    before = deepcopy(matrix.data)
    for method in (:efm, :dfm)
        options = (aperture = [0.01, 0.02], hydraulic_aperture = [0.003, 0.004],
            porosity = [0.3, 0.4])
        model = setup_fractured_reservoir_model(matrix, objects, SinglePhaseSystem();
            method = method, unified = true, subtract_fracture_volume = false,
            fracture_arg = options)
        domain = reservoir_domain(model)
        mesh = physical_representation(domain)
        nm = number_of_cells(mesh.matrix_mesh)
        cells = nm+1:number_of_cells(domain)
        ids = domain[:fracture_id][cells]
        @test domain[:volumes][cells] ≈ tpfv_geometry(mesh.fracture_mesh).volumes.*options.aperture[ids]
        @test domain[:porosity][cells] == options.porosity[ids]
        @test domain[:permeability][1:3, cells] ≈ repeat((options.hydraulic_aperture[ids].^2/12)', 3, 1)
        @test all(iszero, domain[:permeability][4:6, cells])
        parents = [EFM.fracture_parent_cell(mesh.fracture_mesh, mesh.matrix_mesh, c)
            for c in 1:length(cells)]
        @test domain[:region][cells] == domain[:region][parents]
        @test domain[:metadata] == matrix[:metadata]
        domain[:metadata][:values][1] = 10
        @test matrix[:metadata][:values] == [1, 2]
        @test isequal(matrix.data, before)
        explicit = setup_fractured_reservoir_model(matrix, objects, SinglePhaseSystem();
            method = method, unified = true,
            fracture_arg = (options..., permeability = 1e-12))
        @test all(==(1e-12), reservoir_domain(explicit)[:permeability][1:3, cells])
        @test_throws ArgumentError setup_fractured_reservoir_model(matrix, objects,
            SinglePhaseSystem(); method = method, fracture_arg = (hydraulic_aperture = -1,))
        @test_throws ArgumentError setup_fractured_reservoir_model(matrix, objects,
            SinglePhaseSystem(); method = method, fracture_arg = (porosity = 0,))
    end
end

@testset "Face adjacency with isolated cells" begin
    empty = EFM.embedded_cell_faces(zeros(Int, 2, 0), 3)
    @test length(empty) == 3
    @test all(isempty(empty[c]) for c in 1:3)
    connected = EFM.embedded_cell_faces([1 2; 2 3], 4)
    @test [collect(connected[c]) for c in 1:4] == [[1], [1, 2], [2], Int[]]
end

@testset "Small Float32 fracture geometry" begin
    function scaled_fracture_geometry(::Type{T}, scale) where T
        s = T(scale)
        mesh = reservoir_mesh((1, 1, 1), (s, s, s))
        fracture = RectangleFracture(s.*T.([0.4, 0.5, 0.5]), T.([1, 2, 3]), (2s, 2s))
        return embedded_fracture_mesh(mesh, fracture; float_type = T, tolerance = T(1e-5)*s).geometry
    end
    reference = scaled_fracture_geometry(Float64, 1.0)
    scale = 1e-8
    geometry = scaled_fracture_geometry(Float32, scale)
    @test first(geometry.patches).area/scale^2 ≈ first(reference.patches).area rtol = 1e-4
    @test first(geometry.patches).centroid/scale ≈ first(reference.patches).centroid rtol = 1e-4
    @test geometry.matrix_distances/scale ≈ reference.matrix_distances rtol = 1e-4
    matrix = reservoir_domain(reservoir_mesh((1, 1, 1), (scale, scale, scale)))
    fracture = RectangleFracture(Float32(scale).*Float32.([0.5, 0.5, 0.5]), [0, 0, 1],
        (Float32(scale), Float32(scale)); resolution = (3, 4))
    grid = embedded_fracture_mesh(physical_representation(matrix), fracture;
        float_type = Float32, tolerance = Float32(1e-5*scale))
    @test number_of_cells(grid) == 12
    @test sum(tpfv_geometry(grid).volumes)/scale^2 ≈ 1 rtol = 1e-4
    @test EFM.cell_section === Jutul.CutCellMeshes.cell_section
    @test EFM.clip_polygon === Jutul.CutCellMeshes.clip_polygon
    @test EFM.CellSpatialIndex === Jutul.CutCellMeshes.CellSpatialIndex
    @test EFM.FractureFrame === Jutul.CutCellMeshes.PlanarFrame
end
