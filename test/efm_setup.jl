@testset "Parametric embedded geometry" begin
    matrix = reservoir_domain(reservoir_mesh((2, 2, 2), (2.0, 2.0, 2.0)))
    reference = nothing
    for T in (Float64, Float32, BigFloat)
        objects = AbstractEmbeddedFracture[
            RectangleFracture(T.([1, 1, 1]), [1, 0, 0], (2, 2)),
            EllipticalFracture(T.([1, 1, 1]), [0, 1, 0], (0.7, 0.7); segments = 16)]
        @test all(f -> f isa AbstractEmbeddedFracture{T}, objects)
        domain = embedded_fracture_domain(matrix, objects; float_type = T)
        mesh = physical_representation(domain)
        geometry = mesh.geometry
        @test mesh isa EmbeddedFractureMesh{T}
        @test geometry isa EmbeddedFractureGeometry{T}
        @test Jutul.float_type(mesh) == T
        @test eltype(first(geometry.patches).centroid) == T
        @test eltype(geometry.connection_distances) == T
        @test eltype(geometry.matrix_distances) == T
        @test eltype(domain[:volumes]) == T
        @test eltype(domain[:transmissibilities]) == T
        @test all(isfinite, domain[:transmissibilities])
        @test !isempty(geometry.intersections)
        @test eltype(first(geometry.intersections).distances) == T
        combined = combined_fracture_domain(matrix, domain)
        @test combined[:is_fracture] isa Vector{Int}
        @test combined[:is_fracture] == vcat(zeros(Int, number_of_cells(matrix)), ones(Int, number_of_cells(domain)))
        @test physical_representation(combined) isa EmbeddedReservoirMesh{T}
        @test eltype(tpfv_geometry(physical_representation(combined)).volumes) == T
        @test eltype(UnstructuredMesh(physical_representation(combined)).node_points[1]) == T
        plotted = Jutul.triangulate_mesh(physical_representation(combined))
        @test eltype(plotted.points) == T
        if isnothing(reference)
            reference = domain[:volumes]
        else
            @test Float64.(domain[:volumes]) ≈ reference rtol = 1e-5
        end
    end
end

@testset "DFM setup from fracture objects" begin
    matrix = reservoir_domain(reservoir_mesh((2, 2, 2), (2.0, 2.0, 2.0));
        permeability = repeat([2.0, 3.0, 4.0], 1, 8), porosity = collect(0.1:0.01:0.17))
    matrix[:satnum, Cells()] = collect(1:8)
    matrix[:transmissibility_multiplier, Faces()] = collect(2.0:13.0)
    matrix[:net_to_gross, Cells()] = fill(0.4, 8)
    Jutul.set_mesh_entity_tag!(physical_representation(matrix), Faces(), :faults, :test_fault, [1])
    original_vertical = findall(f -> abs(matrix[:normals][3, f]) < 0.5, 1:number_of_faces(matrix))
    Jutul.set_mesh_entity_tag!(physical_representation(matrix), Faces(), :orientation, :vertical, original_vertical)
    before = deepcopy(matrix.data)
    objects = AbstractEmbeddedFracture[
        RectangleFracture([0.4, 0.5, 0.5], [1, 0, 0], (0.7, 0.6); aperture = 0.01, porosity = 0.4),
        EllipticalFracture([1.4, 1.4, 1.4], [1, 1, 1], (0.4, 0.3); segments = 12, aperture = 0.02),
        PolygonFracture([[1.5, 0.1, 0.1], [1.5, 0.9, 0.1], [1.5, 0.9, 0.4],
            [1.5, 0.4, 0.4], [1.5, 0.4, 0.9], [1.5, 0.1, 0.9]]; aperture = 0.03)]
    efm = embedded_fracture_domain(matrix, objects)
    dfm = discrete_fracture_domain(matrix, objects)
    parents, ids = dfm.info[:cell_index], dfm.info[:fracture_id]
    @test dfm.matrix[:satnum] == matrix[:satnum][parents]
    @test dfm.matrix[:porosity] == matrix[:porosity][parents]
    @test dfm.matrix[:permeability] == matrix[:permeability][:, parents]
    @test sum(dfm.matrix[:volumes]) ≈ sum(matrix[:volumes])
    @test all(>(0), dfm.matrix[:volumes])
    @test all(isfinite, dfm.matrix[:cell_centroids])
    mapping = dfm.info[:face_index]
    inherited = findall(>(0), mapping)
    @test dfm.matrix[:transmissibility_multiplier][inherited] == matrix[:transmissibility_multiplier][mapping[inherited]]
    @test all(==(1), dfm.matrix[:transmissibility_multiplier][mapping .== 0])
    @test Jutul.get_mesh_entity_tag(physical_representation(dfm.matrix), Faces(), :faults, :test_fault) == findall(==(1), mapping)
    @test Jutul.mesh_z_is_depth(physical_representation(dfm.matrix)) == Jutul.mesh_z_is_depth(physical_representation(matrix))
    new_faces = findall(==(0), mapping)
    vertical = filter(f -> abs(dfm.matrix[:normals][3, f]) < 1e-10, new_faces)
    horizontal = filter(f -> abs(dfm.matrix[:normals][3, f]) > 1-1e-10, new_faces)
    unit_ntg = deepcopy(dfm.matrix)
    unit_ntg[:net_to_gross, Cells()] .= 1
    base_trans = JutulDarcy.reservoir_transmissibility(unit_ntg)
    trans = JutulDarcy.reservoir_transmissibility(dfm.matrix)
    @test trans[vertical] ≈ 0.4base_trans[vertical]
    @test trans[horizontal] ≈ base_trans[horizontal]
    overridden = deepcopy(matrix)
    overridden[:transmissibility_override, Faces()] = vcat([7.0], fill(NaN, 11))
    override_dfm = discrete_fracture_domain(overridden, objects)
    @test sum(override_dfm.matrix[:transmissibility_override][override_dfm.info[:face_index] .== 1]) ≈ 7.0
    explicit = deepcopy(matrix)
    explicit[:transmissibilities, Faces()] = JutulDarcy.reservoir_transmissibility(matrix)
    explicit_dfm = discrete_fracture_domain(explicit, objects)
    @test sum(explicit_dfm.matrix[:transmissibilities][explicit_dfm.info[:face_index] .== 1]) ≈ explicit[:transmissibilities][1]
    @test all(isfinite, explicit_dfm.matrix[:transmissibilities])
    for id in eachindex(objects)
        @test sum(dfm.fractures[:volumes][ids .== id]) ≈ sum(efm[:volumes][efm[:fracture_id] .== id]) rtol = 1e-7
        @test all(==(objects[id].properties.aperture), dfm.fractures[:aperture][ids .== id])
        @test all(==(objects[id].properties.permeability), dfm.fractures[:permeability][ids .== id])
    end
    fmc = JutulDarcy.FractureMatrixConnection()
    @test all(k -> 2 <= k <= 4, dfm.fractures[:matrix_permeability, fmc])
    @test all(k -> isequal(matrix.data[k], before[k]), keys(before))
    # A small fracture wholly inside one cell must be represented exactly.
    one = reservoir_domain(reservoir_mesh((1, 1, 1)))
    finite = RectangleFracture([0.5, 0.5, 0.5], [1, 0, 0], (0.2, 0.3))
    small = discrete_fracture_domain(one, [finite])
    @test sum(small.fractures[:volumes])/finite.properties.aperture ≈ 0.06
    well = setup_well(one, [1]; name = :Injector, simple_well = true, WI = 1e-10, WIth = 2e-10)
    transferred = EFM.cut_fracture_well(well, small.matrix, small.info[:cell_index])
    @test sum(transferred[:well_index]) ≈ sum(well[:well_index])
    @test sum(transferred[:thermal_well_index]) ≈ sum(well[:thermal_well_index])
    @test sort(physical_representation(transferred).perforations.reservoir) == collect(1:number_of_cells(small.matrix))
    @test transferred[:perforation_centroids] == small.matrix[:cell_centroids]
    for simple in (false, true)
        automatic = setup_well(one, [1]; name = :Automatic, simple_well = simple)
        mapped = EFM.cut_fracture_well(automatic, small.matrix, small.info[:cell_index])
        WI = Jutul.default_parameter_values(automatic, nothing, JutulDarcy.WellIndices(), :WellIndices)
        WIth = Jutul.default_parameter_values(automatic, nothing, JutulDarcy.WellIndicesThermal(), :WellIndicesThermal)
        @test sum(mapped[:well_index]) ≈ sum(WI)
        @test sum(mapped[:thermal_well_index]) ≈ sum(WIth)
        @test all(isnan, automatic[:well_index])
    end
    system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))
    model, parameters = setup_fractured_reservoir_model(one, [finite], system;
        wells = [well], extra_out = true)
    @test haskey(model.models, :Fractures)
    @test haskey(parameters, :Injector)
end

@testset "Same fracture vector in EFM and DFM flow" begin
    matrix = reservoir_domain(reservoir_mesh((4, 4, 3), (4.0, 4.0, 3.0));
        porosity = 0.1, permeability = 0.1si_unit(:darcy))
    objects = [RectangleFracture([2, 2, 1.5], [1, 0, 0], (4, 3); aperture = 1e-4)]
    wells = [setup_well(matrix, [1, 17, 33]; name = :Injector, simple_well = true),
        setup_well(matrix, [16, 32, 48]; name = :Producer, simple_well = true)]
    system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()); reference_densities = (1000.0, 1000.0))
    efm = setup_fractured_reservoir_model(matrix, Any[objects...], system;
        method = :efm, wells = wells)
    dfm = setup_fractured_reservoir_model(matrix, Any[objects...], system;
        method = :dfm, wells = wells)
    control = Dict(:Injector => InjectorControl(TotalRateTarget(1e-5), [1.0, 0.0]; density = 1000.0),
        :Producer => ProducerControl(BottomHolePressureTarget(50si_unit(:bar))))
    results = []
    for model in (efm, dfm)
        state0 = setup_reservoir_state(model; Pressure = 50si_unit(:bar), Saturations = [0.0, 1.0])
        forces = setup_reservoir_forces(model; control = control)
        case = JutulCase(model, fill(200.0, 3), forces; state0 = state0)
        solver = GenericKrylov(preconditioner = ILUZeroPreconditioner(), relative_tolerance = 1e-10)
        push!(results, simulate_reservoir(case; info_level = -1, linear_solver = solver,
            tol_cnv = 1e-10, timesteps = :none))
    end
    @test length(results[1].result.states) == length(results[2].result.states) == 3
    for (e, d) in zip(results[1].result.states, results[2].result.states)
        @test e[:Reservoir][:Pressure] ≈ d[:Reservoir][:Pressure] rtol = 1e-9
        @test e[:Reservoir][:Saturations] ≈ d[:Reservoir][:Saturations] atol = 1e-9
    end
end
