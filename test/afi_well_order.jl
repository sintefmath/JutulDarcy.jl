using Dates, GeoEnergyIO, Jutul, JutulDarcy, Test

@testset "AFI perforation order and masks" begin
    ix = GeoEnergyIO.IXParser
    afi = JutulDarcy.AFISetup
    function records(texts...)
        parsed = ix.parse_ix_record.(texts)
        units = ix.get_unit_system_pair(:si, :si; ix_dict = ix.conversion_ix_dict())
        return ix.convert_ix_records(collect(parsed), "test", units)
    end
    function connection_update(fields...)
        return (keyword = "WellDef", value = Dict(
            "WellName" => "P", "WellToCellConnections" => Dict(fields...)))
    end
    date = DateTime(2026, 1, 1)
    days = date .+ Day.(0:6)
    ix_steps = Jutul.OrderedDict(
        days[1] => [],
        days[2] => vcat(records("""
            WellDef "P" {
                WellToCellConnections [
                    Cell Status PiMultiplier
                    100 OPEN 0.4
                    105 OPEN 0.8
                ]
            }
            """), [connection_update("Status" => [missing, missing, ix.IX_CLOSED],
                "PiMultiplier" => [missing, missing, 0.6])]),
        days[3] => vcat([connection_update("Status" => [missing, missing, ix.IX_OPEN])], records("""
            WellDef "P" {
                WellToCellConnections [
                    Cell Status PiMultiplier WellBoreRadius
                    102 OPEN 0.25 0.23
                ]
            }
            """)),
        days[4] => [],
        days[5] => [],
        days[6] => [],
        days[7] => [],
    )
    fm_steps = Jutul.OrderedDict(
        days[1] => records("""
            Well "P" {
                BottomHoleRefDepth=50.0
                Constraints[BOTTOM_HOLE_PRESSURE]=1.0e7
            }
            """),
        days[2] => [],
        days[3] => [],
        days[4] => records("""Well "P" { Status=ALL_COMPLETIONS_SHUTIN }"""),
        days[5] => records("""Well "P" { Constraints[BOTTOM_HOLE_PRESSURE]=2.0e7 }"""),
        days[6] => records("""Well "P" { Status=OPEN }"""),
        days[7] => [],
    )
    ix_model = records("""WellDef "Q" { Undefined=TRUE }""", """
        WellDef "P" {
            WellToCellConnections [
                Cell Completion PenetrationDirection Status PiMultiplier Transmissibility
                105 "deep" K CLOSED 0.7 3.0e-12
                101 "inactive" K OPEN 0.9 9.0e-12
                100 "shallow" K OPEN 0.5 1.0e-12
            ]
        }
        """, """
        WellDef "Q" {
            WellToCellConnections [
                Cell Completion PenetrationDirection Status
                104 "A" I CLOSED
                104 "B" I OPEN
                100 "C" I OPEN
            ]
        }
        """)
    pushfirst!(ix_model, (keyword = "StructuredInfo", value = Dict("FirstCellId" => 100)))
    input = ix.AFIInputFile(Dict(
        "IX" => Dict("MODEL_DEFINITION" => ix_model, "STEPS" => ix_steps),
        "FM" => Dict("MODEL_DEFINITION" => [], "STEPS" => fm_steps),
    ))
    mesh = extract_submesh(CartesianMesh((2, 1, 3), (20.0, 10.0, 30.0)), [1, 3, 4, 5, 6])
    domain = reservoir_domain(mesh; permeability = 1.0e-13, porosity = 0.2)
    perfs = JutulDarcy.Perforations()
    sys = SinglePhaseSystem(AqueousPhase())

    for (options, permutation, q_order) in (
            ((perf_sort = Dict("P" => :depth, "Q" => :track),), [1, 2, 3], [3, 1, 2]),
            (NamedTuple(), [3, 1, 2], [1, 2, 3]))
        wells = afi.setup_wells(input, domain; options...)
        q, p = wells
        @test physical_representation(p).perforations.reservoir == [1, 2, 5][permutation]
        @test p[:original_perforation_indices, perfs] == [3, 4, 1][permutation]
        @test p[:perforation_connection_ids, perfs] == [(100, "shallow"), (102, nothing), (105, "deep")][permutation]
        @test p[:perforation_radius, perfs] == [0.1, 0.23, 0.1][permutation]
        @test p[:cell_centroids, Cells()][3, 1] == 50.0
        @test isequal(p[:well_index, perfs], [1.0e-12, NaN, 3.0e-12][permutation])
        @test q[:perforation_names, perfs] == ["A", "B", "C"][q_order]
        model = setup_reservoir_model(domain, sys; wells = wells)
        @test collect(keys(JutulDarcy.get_model_wells(model))) == [:Q, :P]
        @test_throws DimensionMismatch setup_forces(model[:P], mask = PerforationMask([1.0, 0.0]))
        flow = only(term.cross_term for term in model.cross_terms if term.source == :P && term.target == :Reservoir)
        @test flow.reservoir_cells == [1, 2, 5][permutation]
        _, forces = afi.setup_afi_schedule(input, model)
        @test [f[:P].mask.values for f in forces] == [
            [0.5, 0.0, 0.0][permutation],
            [0.0, 0.0, 0.8][permutation],
            [0.6, 0.25, 0.8][permutation],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.6, 0.25, 0.8][permutation],
        ]
        @test all(f -> f[:Q].mask.values == [0.0, 1.0, 1.0][q_order], forces)
        @test all(f -> f[:Facility].control[:P] isa JutulDarcy.DisabledControl, forces[4:5])
    end
end
