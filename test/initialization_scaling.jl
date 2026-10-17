using Test, Jutul, GeoEnergyIO, JutulDarcy
import JutulDarcy as JD

@testset "SWATINIT preserves saturation closure" begin
    data = parse_data_file(joinpath(@__DIR__,"data","swatinit.DATA");silent=true)
    case = setup_case_from_data_file(data)
    S = case.state0[:Reservoir][:Saturations]
    @test all(isapprox.(vec(sum(S;dims=1)),1.0))
    @test S[:,11] ≈ [0.3,0.7]
    @test S[1,3] ≈ 0.2
end

@testset "Three-phase scaling uses table endpoints" begin
    krw = JD.PhaseRelativePermeability([0.,0.2,1.],[0.,0.,1.])
    kro = JD.PhaseRelativePermeability([0.,1.],[0.,1.])
    sc(k) = reshape(collect(JD.get_kr_scalers(k)),4,1)
    tables = (krw,kro,kro,kro)
    scaled = JD.get_three_phase_scaled_relperms(JD.ThreePointKrScale(),tables...,0.,sc(krw),sc(kro),sc(kro),sc(kro),1)
    for (base,curve) in zip(tables,scaled)
        @test maximum(abs(curve(s)-base(s)) for s in 0:0.01:1) <= 2e-4
    end
    # Changing the cell's connate water shifts the scaled saturation axis.
    # Table-side residual endpoints must stay tied to the unchanged tables.
    sw = sc(krw)
    sw[1] = 0.1
    _,_,og,g = JD.get_three_phase_scaled_relperms(JD.ThreePointKrScale(),tables...,0.1,sw,sc(kro),sc(kro),sc(kro),1)
    @test og(0.9) ≈ kro(1.0)
    @test g(0.9) ≈ kro(1.0)
    @test og(0.45) ≈ kro(0.5) atol=2e-4
end
