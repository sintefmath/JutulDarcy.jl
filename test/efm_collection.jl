struct ShiftedFractureVector{T, V <: AbstractVector{T}} <: AbstractVector{T}
    values::V
end
Base.size(v::ShiftedFractureVector) = size(v.values)
Base.axes(v::ShiftedFractureVector) = (0:length(v.values)-1,)
Base.getindex(v::ShiftedFractureVector, i::Int) = v.values[i+1]

@testset "Fracture collections and general vectors" begin
    r1 = RectangleFracture([0.5, 1, 1], [1, 0, 0], (2, 2); aperture = 0.01)
    e = EllipticalFracture([1, 1, 1], [0, 1, 0], (0.7, 0.8); segments = 16)
    p = PolygonFracture([[1.5, 0, 0], [1.5, 1, 0], [1.5, 1, 1], [1.5, 0, 1]])
    r2 = RectangleFracture([1.8, 1, 1], [1, 0, 0], (2, 2); aperture = 0.02)
    r32 = RectangleFracture(Float32.([0.2, 1, 1]), [1, 0, 0], (2, 2))
    raw = Any[r1, e, p, r2, r32]
    grouped = FractureCollection(raw)
    @test length(grouped) == 5
    @test length(grouped.groups) == 4
    @test all(g -> isconcretetype(eltype(g)), grouped.groups)
    @test grouped.indices == ([1, 4], [2], [3], [5])
    @test all(i -> grouped[i] === raw[i], eachindex(raw))
    @test FractureCollection(grouped) === grouped
    @test collect(FractureCollection(ShiftedFractureVector(raw))) == raw
    @test collect(FractureCollection(@view raw[2:4])) == raw[2:4]
    @test length(FractureCollection(r1)) == 1
    @test isempty(FractureCollection(Any[]))
    @test_throws ArgumentError FractureCollection(Any[r1, 1])
    @test_throws ArgumentError FractureCollection([1, 2])
    @test_throws BoundsError grouped[0]
    matrix = reservoir_domain(reservoir_mesh((2, 2, 2), (2.0, 2.0, 2.0)))
    domain = embedded_fracture_domain(matrix, grouped; aperture = collect(0.01:0.01:0.05))
    mesh = physical_representation(domain)
    ids = domain[:fracture_id]
    @test mesh.fractures isa FractureCollection
    @test issorted(ids)
    @test domain[:aperture] ≈ collect(0.01:0.01:0.05)[ids]
    @test all(>(0), domain[:volumes])
    # Regrouping must retain each object's geometry and all shared-face paths.
    for id in eachindex(raw)
        individual = embedded_fracture_domain(matrix, raw[id])
        @test domain[:cell_centroids][:, ids .== id] ≈ individual[:cell_centroids]
        @test tpfv_geometry(mesh).volumes[ids .== id] ≈ tpfv_geometry(physical_representation(individual)).volumes
    end
    for value in (grouped, mesh, mesh.geometry, physical_representation(combined_fracture_domain(matrix, domain)))
        output = sprint(show, MIME"text/plain"(), value)
        @test occursin(string(nameof(typeof(value))), output)
        @test ncodeunits(output) < 1000
    end
    @test occursin("RectangleFracture{Float32}", sprint(show, MIME"text/plain"(), grouped))
    @test occursin("4 vertices", sprint(show, p))
end

@testset "Unified EFM and DFM setup" begin
    matrix = reservoir_domain(reservoir_mesh((2, 2, 2), (2.0, 2.0, 2.0)))
    fracture = RectangleFracture([1, 1, 1], [1, 0, 0], (2, 2); aperture = 0.01)
    system = ImmiscibleSystem((AqueousPhase(), LiquidPhase()))
    original = sum(matrix[:volumes])
    for objects in ([fracture], Any[fracture], AbstractEmbeddedFracture[fracture],
            @view(Any[fracture, fracture][1:1]), ShiftedFractureVector(Any[fracture]),
            FractureCollection(Any[fracture]))
        for method in (:efm, :dfm)
            model, parameters = setup_fractured_reservoir_model(matrix, objects, system;
                method = method, extra_out = true)
            @test haskey(model.models, :Fractures)
            @test haskey(parameters, :Fractures)
            @test sum(reservoir_domain(model)[:volumes]) + sum(model[:Fractures].data_domain[:volumes]) ≈ original
            @test physical_representation(model[:Fractures].data_domain) isa Jutul.EmbeddedMesh || method == :efm
            displayed = sprint(show, MIME"text/plain"(), physical_representation(model[:Fractures].data_domain))
            @test occursin("cells", displayed)
            @test ncodeunits(displayed) < 1000
        end
    end
    single = setup_fractured_reservoir_model(matrix, Any[fracture], system;
        method = :efm, unified = true)
    @test physical_representation(reservoir_domain(single)) isa EmbeddedReservoirMesh
    @test !haskey(single.models, :Fractures)
    domain = reservoir_domain(single)
    @test domain[:is_fracture, Cells()] isa Vector{Int}
    @test domain[:is_fracture] == vcat(zeros(Int, number_of_cells(matrix)), ones(Int, 4))
    embedded = embedded_fracture_domain(matrix, fracture)
    preprocessed = setup_fractured_reservoir_model(matrix, embedded, system; unified = true)
    @test reservoir_domain(preprocessed)[:is_fracture] == domain[:is_fracture]
    @test !haskey(preprocessed.models, :Fractures)
    for method in (:efm, :dfm)
        model = setup_fractured_reservoir_model(matrix, Any[fracture], system;
            method = method, unified = false, subtract_fracture_volume = false)
        @test haskey(model.models, :Fractures)
        @test sum(reservoir_domain(model)[:volumes]) ≈ original
        @test_throws ArgumentError setup_fractured_reservoir_model(matrix, Any[fracture, 1], system; method = method)
        @test_throws ArgumentError setup_fractured_reservoir_model(matrix, Any[], system; method = method)
    end
    @test_throws ArgumentError setup_fractured_reservoir_model(matrix, [fracture], system; method = :bad)
    dfm_single = setup_fractured_reservoir_model(matrix, Any[fracture], system;
        method = :dfm, unified = true)
    @test !haskey(dfm_single.models, :Fractures)
    @test reservoir_domain(dfm_single)[:is_fracture] == domain[:is_fracture]
    @test_throws TypeError setup_fractured_reservoir_model(matrix, [fracture], system;
        method = :efm, unified = :bad)
    efm_extended = setup_fractured_reservoir_model(matrix, [fracture], system;
        method = :efm, extend_wells = true)
    @test haskey(efm_extended.models, :Fractures)
    @test sum(matrix[:volumes]) == original
end
