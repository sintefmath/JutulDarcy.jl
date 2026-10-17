using Test, Jutul, GeoEnergyIO, JutulDarcy
import JutulDarcy as JD

@testset "Schedule completion interpretation" begin
    sys = ImmiscibleSystem((AqueousPhase(),LiquidPhase()))
    @test_throws "THP control not supported" JD.producer_control(sys,"OPEN","THP",1.,0.,0.,1.,1e7)
    dom = reservoir_domain(CartesianMesh((1,1,1),(100.,100.,10.));permeability=1e-13,porosity=0.2)
    for steps in ([],[Dict{String,Any}()])
        schedule = Dict("WELSPECS"=>Dict(),"STEPS"=>steps)
        @test_throws "SCHEDULE contains no report steps" JD.parse_schedule(dom,Dict(),Dict(),schedule,sys)
    end
    perf = (WI=-1.,diameter=0.2,Kh=NaN,skin=0.,dir="Z",mul=1.,ctrl=1,open=true,r0=5.)
    v = Dict((1,1,1)=>perf)
    calc(wspec) = JD.compdat_to_connection_factors(dom,wspec,v,1;sort=false)[2][1]
    @test calc((drainage_radius=-1.,)) ≈ 2π*1e-13*10/log(5/0.1)
    @test calc((drainage_radius=100.,)) ≈ calc((drainage_radius=-1.,))
    v[(1,1,1)] = merge(perf,(r0=-1.,))
    @test calc((drainage_radius=100.,)) ≈ JD.compute_peaceman_index(physical_representation(dom),1e-13,0.1,1)
    v[(1,1,1)] = merge(perf,(WI=1e-10,))
    @test calc((drainage_radius=-1.,)) == 1e-10

    data = parse_data_file(GeoEnergyIO.test_input_file_path("SPE1","SPE1.DATA");silent=true)
    base = deepcopy(first(data["SCHEDULE"]["STEPS"]))
    base["TSTEP"] = [86400.]
    base["WPIMULT"] = [["PROD",3.,-1,-1,-1,-1,-1]]
    later = Jutul.OrderedDict{String,Any}("WCONPROD"=>base["WCONPROD"],"TSTEP"=>[86400.])
    again = deepcopy(later)
    again["WPIMULT"] = [["PROD",2.,-1,-1,-1,-1,-1]]
    reset = deepcopy(later)
    reset["COMPDAT"] = deepcopy(base["COMPDAT"])
    data["SCHEDULE"]["STEPS"] = [base,deepcopy(later),again,deepcopy(later),reset]
    case = setup_case_from_data_file(data)
    multiplier(f) = isnothing(f[:PROD].mask) ? 1.0 : only(f[:PROD].mask.values)
    @test multiplier.(case.forces) ≈ [3.,3.,6.,6.,1.]
end
