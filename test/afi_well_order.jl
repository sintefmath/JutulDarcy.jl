using Dates, GeoEnergyIO, Jutul, JutulDarcy, Test

@testset "AFI depth order and perforation masks" begin
    ix = GeoEnergyIO.IXParser
    afi = JutulDarcy.AFISetup
    open, closed = ix.IX_OPEN, ix.IX_CLOSED
    welldef(fields...) = (keyword = "WellDef", value = Dict(
        "WellName" => "P", "WellToCellConnections" => Dict{String, Any}(fields...)))
    initial = welldef(
        "Cell" => [105, 101, 100, 100],
        "Completion" => ["B", "inactive", "A", "B"],
        "TrueVerticalDepth" => [25.0, 2.0, 5.0, 4.0],
        "Status" => [open, open, closed, open],
        "PiMultiplier" => [0.8, 0.9, 0.7, 0.5],
    )
    reordered = welldef(
        "Cell" => [100, 105, 100],
        "Completion" => ["A", "B", "B"],
        "Status" => [open, closed, open],
        "PiMultiplier" => [0.1, 0.2, 0.3],
        "WellBoreRadius" => [0.16, 0.18, 0.2],
    )
    indexed = welldef("Status" => [missing, missing, closed])
    legacy_indexed = welldef("Status" => [open, open, open], "PiMultiplier" => [1.0, 1.0, 1.0])
    control = (keyword = "Well", value = Dict(
        "name" => "P", "Constraints" => ix.IXAssignmentRecord("Constraints", "BOTTOM_HOLE_PRESSURE", 1.0e7)))
    dates = DateTime(2026, 1, 1) .+ Hour.(0:3)
    ix_steps = Jutul.OrderedDict(zip(dates, ([], [reordered], [indexed, legacy_indexed], [])))
    fm_steps = Jutul.OrderedDict(zip(dates, ([control], [], [], [])))
    offset = (keyword = "StructuredInfo", value = Dict("FirstCellId" => 100))
    input = ix.AFIInputFile(Dict(
        "IX" => Dict("MODEL_DEFINITION" => [offset, initial], "STEPS" => ix_steps),
        "FM" => Dict("MODEL_DEFINITION" => [], "STEPS" => fm_steps),
    ))

    mesh = extract_submesh(CartesianMesh((2, 1, 3), (20.0, 10.0, 30.0)), [1, 3, 4, 5, 6])
    domain = reservoir_domain(mesh; permeability = 1.0e-13, porosity = 0.2)
    wells = afi.setup_wells(input, domain)
    well = only(wells)
    perfs = JutulDarcy.Perforations()
    @test physical_representation(well).perforations.reservoir == [1, 1, 5]
    @test well[:original_perforation_indices, perfs] == [4, 3, 1]
    @test well[:perforation_connection_ids, perfs] == [(100, "B"), (100, "A"), (105, "B")]
    @test well[:perforation_centroids, perfs][3, :] == [4.0, 5.0, 25.0]
    @test well[:perforation_radius, perfs] == [0.2, 0.16, 0.18]

    sys = SinglePhaseSystem(AqueousPhase())
    model = setup_reservoir_model(domain, sys; wells = wells, backend = :csc, block_backend = false)
    dt, forces = afi.setup_afi_schedule(input, model)
    @test [f[:P].mask.values for f in forces] == [
        [0.5, 0.0, 0.8], [0.3, 0.1, 0.0], [0.3, 0.0, 0.0]
    ]
    state0 = setup_reservoir_state(model; Pressure = 1.2e7)
    case = JutulCase(model, dt, forces; state0 = state0)
    result = simulate_reservoir(case; info_level = -1, linear_solver = nothing)
    @test length(result.states) == length(dt)
    @test all(<(0.0), result.wells[:P][:mass_rate])
end
