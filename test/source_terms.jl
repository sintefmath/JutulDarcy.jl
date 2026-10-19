using Test, Jutul, JutulDarcy

@testset "Source types and concrete source arrays" begin
    domain = get_1d_reservoir(2; L=1.0,poro=0.1,area=1.0)
    sys = ImmiscibleSystem((LiquidPhase(),VaporPhase());reference_densities=[800.0,2.0])
    model = SimulationModel(domain,sys)
    state = (RelativePermeabilities=[1.0 1.0;0.25 0.25],
        PhaseViscosities=ones(2,2), PhaseMassDensities=[1000.0 1000.0;10.0 10.0])
    eq = model.equations[:mass_conservation]
    storage = (state=state,)
    for (type,density) in ((MassSource,[1.0,1.0]),(StandardVolumeSource,[800.0,2.0]),(VolumeSource,[1000.0,10.0]))
        src = SourceTerm(1,0.1;fractional_flow=[0.25,0.75],type=type)
        acc = zeros(2,2)
        Jutul.apply_forces_to_equation!(acc,storage,model,eq,nothing,[src],0.0)
        @test acc[:,1] ≈ -0.1 .* [0.25,0.75] .* density
        @test iszero(acc[:,2])
        @test only(Jutul.subforce([src],model)).type == type
        src = SourceTerm(1,-0.1;type=type)
        fill!(acc,0)
        Jutul.apply_forces_to_equation!(acc,storage,model,eq,nothing,[src],0.0)
        @test acc[:,1] ≈ 0.1 .* [0.8,0.2] .* density
    end
    src = SourceTerm(1,0.1;fractional_flow=[1.0,0.0])
    @test setup_forces(model;sources=src).sources == [src]
    for sources in (SourceTerm[src], SourceTerm[], Any[src])
        @test_throws ArgumentError setup_forces(model;sources=sources)
    end
    @test_throws ArgumentError Jutul.apply_forces_to_equation!(zeros(2,2),storage,model,eq,nothing,SourceTerm[src],0.0)
end
