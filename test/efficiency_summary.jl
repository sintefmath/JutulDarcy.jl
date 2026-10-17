using Test, Jutul, GeoEnergyIO, JutulDarcy
import JutulDarcy as JD

@testset "Efficiency-weighted summaries" begin
    dom = reservoir_domain(CartesianMesh((1,1,1),(100.,100.,10.));permeability=1e-13,porosity=0.2)
    sys = ImmiscibleSystem((AqueousPhase(),LiquidPhase()))
    model = MultiModel((Reservoir=SimulationModel(dom,sys),))
    wi = Dict{Symbol,AbstractVector}(:wrat=>[2.,2.],:orat=>[0.,0.],:bhp=>[1e7,1e7],:factor=>[0.8,0.5])
    wp = Dict{Symbol,AbstractVector}(:wrat=>[-1.,-1.],:orat=>[-3.,-3.],:bhp=>[1e7,1e7],:factor=>[0.7,0.9])
    wells = JD.WellResults([10.,20.],Dict(:I=>wi,:P=>wp),nothing)
    summary = JD.summary_result(model,wells,missing,:si)
    vals = summary["VALUES"]
    @test vals["FIELD"]["FWIR"] ≈ [1.6,1.]
    @test vals["FIELD"]["FWIT"] ≈ [16.,26.]
    @test vals["FIELD"]["FOPR"] ≈ [2.1,2.7]
    @test vals["FIELD"]["FWCT"] ≈ [0.25,0.25]
    @test vals["WELLS"]["I"]["WWIR"] == [2.,2.]
    @test vals["WELLS"]["I"]["WWIT"] ≈ [16.,26.]
    @test vals["WELLS"]["P"]["WOPT"] ≈ [21.,48.]
    @test vals["WELLS"]["I"]["WEFAC"] == [0.8,0.5]
    group = JD.reservoir_measurables(model,wells;type=:group,wells=[:I])
    @test group[:gwir].values ≈ [1.6,1.]
    recomputed = JD.expand_summary(summary;recompute=true)["VALUES"]
    @test recomputed["FIELD"]["FWIR"] ≈ [1.6,1.]
    @test recomputed["WELLS"]["I"]["WWIT"] ≈ [16.,26.]
    delete!(wi,:factor)
    @test JD.reservoir_measurables(model,wells;type=:well,wells=:I)[:wwit].values ≈ [20.,40.]

    # Factors are extracted from the same step's well control.
    case = setup_case_from_data_file(GeoEnergyIO.test_input_file_path("SPE1","SPE1.DATA"))
    f1,f2 = deepcopy(case.forces[1]),deepcopy(case.forces[1])
    for (force,factor) in ((f1,0.8),(f2,0.5))
        old = force[:Facility].control[:PROD]
        force[:Facility].control[:PROD] = ProducerControl(old.target;factor=factor)
    end
    @test JD.well_output(case.model,[case.state0,case.state0],:PROD,[f1,f2],:factor) == [0.8,0.5]
end
