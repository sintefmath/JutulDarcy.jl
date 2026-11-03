using Jutul, JutulDarcy, Test

function group_summary_test_case()
    domain = reservoir_domain(CartesianMesh((5, 1, 1), (500.0, 100.0, 20.0));
        permeability = 1e-13, porosity = 0.2)
    names = [:I1, :P1, :I2, :P2, :Other]
    wells = []
    for i in eachindex(names)
        push!(wells, setup_vertical_well(domain, i, 1; name = names[i], simple_well = true))
    end
    system = ImmiscibleSystem((AqueousPhase(), LiquidPhase(), VaporPhase()))
    groups = Dict(:Injection => [:I1, :I2], :Production => [:P1, :P2])
    model = setup_reservoir_model(domain, system; wells, well_groups = groups)
    data = Dict(
        :I1 => Dict(:wrat => [2.0, 2.0], :orat => [1.0, 1.0], :grat => [3.0, 3.0],
            :bhp => [1e7, 1e7], :factor => [0.8, 0.5]),
        :I2 => Dict(:wrat => [4.0, 4.0], :orat => [2.0, 2.0], :grat => [6.0, 6.0],
            :bhp => [1e7, 1e7], :factor => [0.25, 0.75]),
        :P1 => Dict(:wrat => [-1.0, -1.0], :orat => [-3.0, -3.0], :grat => [-8.0, -8.0],
            :bhp => [1e7, 1e7], :factor => [0.5, 0.25]),
        :P2 => Dict(:wrat => [-6.0, -6.0], :orat => [-2.0, -2.0], :grat => [-4.0, -4.0],
            :bhp => [1e7, 1e7], :factor => [0.75, 0.5]),
        :Other => Dict(:wrat => [-2.0, -2.0], :orat => [-1.0, -1.0], :grat => [0.0, 0.0],
            :bhp => [1e7, 1e7])
    )
    return model, JutulDarcy.WellResults([10.0, 30.0], data)
end

@testset "Named group summaries" begin
    model, wells = group_summary_test_case()
    summary = JutulDarcy.summary_result(model, wells, missing, :si)
    groups = summary["VALUES"]["GROUP"]
    @test Set(keys(groups)) == Set(["Injection", "Production"])
    injection = groups["Injection"]
    production = groups["Production"]
    @test injection["GWIR"] ≈ [2.6, 4.0]
    @test injection["GOIR"] ≈ [1.3, 2.0]
    @test injection["GGIR"] ≈ [3.9, 6.0]
    @test injection["GLIR"] ≈ [3.9, 6.0]
    @test injection["GWIT"] ≈ [26.0, 106.0]
    @test injection["GLIT"] ≈ [39.0, 159.0]
    @test production["GWPR"] ≈ [5.0, 3.25]
    @test production["GOPR"] ≈ [3.0, 1.75]
    @test production["GGPR"] ≈ [7.0, 4.0]
    @test production["GLPR"] ≈ [8.0, 5.0]
    @test production["GWPT"] ≈ [50.0, 115.0]
    @test production["GOPT"] ≈ [30.0, 65.0]
    @test production["GWCT"] ≈ [0.625, 0.65]
    @test production["GGOR"] ≈ [7/3, 4/1.75]
    @test all(iszero, injection["GWPR"])
    @test all(iszero, production["GWIR"])
    @test summary["VALUES"]["FIELD"]["FWPR"] ≈ production["GWPR"] .+ 2.0

    metric = JutulDarcy.summary_result(model, wells, missing, :metric)
    @test metric["VALUES"]["GROUP"]["Injection"]["GWIR"] ≈ injection["GWIR"]*si_unit(:day)
    @test metric["VALUES"]["GROUP"]["Production"]["GWPT"] ≈ production["GWPT"]
    @test metric["VALUES"]["GROUP"]["Production"]["GWCT"] ≈ production["GWCT"]
    @test metric["VALUES"]["WELLS"]["P1"]["WOPT"] ≈ [15.0, 30.0]
    @test metric["VALUES"]["FIELD"]["FWPT"] ≈ [70.0, 175.0]
    group_only = JutulDarcy.summary_result(model, wells, missing, :si;
        wells = false, field = false, groups = true)
    @test group_only["VALUES"]["GROUP"] == groups
    @test group_only["TIME"].seconds == wells.time
    @test !haskey(group_only["VALUES"], "WELLS")
    @test !haskey(JutulDarcy.summary_result(model, wells; groups = false)["VALUES"], "GROUP")

    # Group vectors from flat summary files retain the same nested representation.
    flat = Dict("YEARS" => wells.time/si_unit(:year),
        "GWIR" => Dict("Injection" => injection["GWIR"]),
        "WWIR" => Dict("I1" => wells.wells[:I1][:wrat]))
    converted = JutulDarcy.convert_summary(flat; unit_system = "SI")
    @test converted["VALUES"]["GROUP"]["Injection"]["GWIR"] == injection["GWIR"]
    @test converted["VALUES"]["WELLS"]["I1"]["WWIR"] == [2.0, 2.0]
end
