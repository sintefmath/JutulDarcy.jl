using Test, JutulDarcy, Jutul

function test_kr_bounds(relperm, n)
    for ph in 1:n
        S = zeros(n, 1)
        S[ph] = 1.0

        kr = similar(S)
        JutulDarcy.update_kr!(kr, relperm, nothing, S, entity_eachindex(kr))
        for i in 1:n
            if i == ph
                @test kr[i] > 0
            else
                @test kr[i] == 0
            end
        end
    end
end

function kr_test_sat(n)
    if n == 2
        S = zeros(2, 3)
        S[1, 1] = 1.0
        S[1, 2] = 0.5
        S[1, 3] = 0.1
        for i = 1:3
            S[2, i] = 1 - S[1, i]
        end
    else
        error("Bad input $n")
    end
    return S
end

function test_brooks_corey_kr()
    bc = BrooksCoreyRelativePermeabilities(2, [2.0, 3.0], [0.2, 0.3], [0.9, 1.0])
    S = kr_test_sat(2)
    kr = similar(S)
    @test JutulDarcy.update_kr!(kr, bc, nothing, S, entity_eachindex(kr)) ≈ [0.9 0.324 0.0; 0.0 0.064 1.0]
    test_kr_bounds(bc, 2)
end

function test_standard_kr()
    S = kr_test_sat(2)
    kr = similar(S)

    kr_1 = S -> S^2
    kr_2 = S -> S
    rel_single = JutulDarcy.RelativePermeabilities((kr_1, kr_2))
    JutulDarcy.update_kr!(kr, rel_single, nothing, S, entity_eachindex(kr))
    for i in axes(kr, 2)
        @test kr[1, i] == kr_1(S[1, i])
        @test kr[2, i] == kr_2(S[2, i])
    end

    # Quadratic in first region, linear in second
    rel_regs = JutulDarcy.RelativePermeabilities(((kr_1, kr_2), (kr_1, kr_2)), regions = [1, 2])
    S = repeat([0.5], 2, 2)
    kr = similar(S)
    JutulDarcy.update_kr!(kr, rel_regs, nothing, S, entity_eachindex(kr))
    @test kr[1, 1] == kr[2, 1] ≈ 0.25
    @test kr[1, 2] == kr[2, 2] ≈ 0.5
end

function test_rel_perm_wrapper()
    s = [0.1, 0.15, 0.2, 0.8, 1.0]
    kr = [0.0, 0.0, 0.4, 0.9, 0.9]
    relperm = PhaseRelativePermeability(s, kr)
    @testset "Detection of points" begin
        @test relperm.connate == 0.1
        @test relperm.k_max == 0.9
        @test relperm.s_max == 0.8
        @test relperm.critical == 0.15
    end
    kro = @. (1 - s)^2
    swof = hcat(s, kr, kro)
    krw, krow = table_to_relperm(swof)
    @testset "SWOF-like table" begin
        @test krw.(s) ≈ swof[:, 2]
        @test krow.(1 .- s) ≈ swof[:, 3]
    end
end

@testset "RelativePermeabilities" begin
    @testset "BrooksCorey" begin
        test_brooks_corey_kr()
    end
    @testset "Standard kr" begin
        test_standard_kr()
    end
    @testset "Wrapper" begin
        test_rel_perm_wrapper()
    end
end

##
using Test, JutulDarcy, Jutul
function solve_1d_kr(phases;
        nc = 3,
        scaling = JutulDarcy.NoKrScale(),
        non_wetting_hysteresis = JutulDarcy.NoHysteresis(),
        wetting_hysteresis = JutulDarcy.NoHysteresis(),
        time = 1.0,
        nstep = nc
    )
    T = time
    tstep = repeat([T/nstep], 2*nstep)
    domain = get_1d_reservoir(nc)
    nc = number_of_cells(domain)
    timesteps = tstep*3600*24
    bar = 1e5
    p0 = 100*bar

    s = collect(range(0, 1, 100))
    drainage = PhaseRelativePermeability(s, s)
    imbibition = PhaseRelativePermeability(s, s.^3)
    should_use_hysteresis = !(wetting_hysteresis == JutulDarcy.NoHysteresis() && non_wetting_hysteresis == JutulDarcy.NoHysteresis())
    if should_use_hysteresis
        relperm = (drainage, imbibition)
    else
        relperm = drainage
    end

    krw = nothing
    krg = nothing
    krog = nothing
    krow = nothing
    all_scalers = [:RelPermScalingW, :RelPermScalingOW, :RelPermScalingOG, :RelPermScalingG]
    if phases == :wog
        mu = [1e-3, 5e-3, 1e-4]
        f1 = [0.5, 0.5, 0.0]
        f2 = [0.0, 0.0, 1.0]
        s0 = [0.2, 0.0, 0.8]
        mph = (AqueousPhase(), LiquidPhase(), VaporPhase())
        krow = krog = krw = krg = relperm
        scalers = all_scalers
    else
        mu = [1e-3, 5e-3]
        f1 = [1.0, 0.0]
        f2 = [0.0, 1.0]
        s0 = [0.0, 1.0]
        if phases == :wo
            krw = relperm
            krow = relperm
            mph = (AqueousPhase(), LiquidPhase())
            scalers = [:RelPermScalingW, :RelPermScalingOW]
        elseif phases == :wg
            krw = relperm
            krg = relperm
            mph = (AqueousPhase(), VaporPhase())
            scalers = [:RelPermScalingW, :RelPermScalingG]
        elseif phases == :og
            krog = relperm
            krg = relperm
            mph = (LiquidPhase(), VaporPhase())
            scalers = [:RelPermScalingOG, :RelPermScalingG]
        else
            error("Bad phases $phases")
        end
    end
    sys = ImmiscibleSystem(mph)
    model = SimulationModel(domain, sys)
    kr = JutulDarcy.ReservoirRelativePermeabilities(
        w = krw,
        g = krg,
        ow = krow,
        og = krog,
        scaling = scaling,
        hysteresis_w = wetting_hysteresis,
        hysteresis_ow = non_wetting_hysteresis,
        hysteresis_og = non_wetting_hysteresis,
        hysteresis_g = non_wetting_hysteresis
    )
    replace_variables!(model, RelativePermeabilities = kr)
    JutulDarcy.add_relperm_parameters!(model)
    tot_time = sum(timesteps)
    deps = Jutul.get_dependencies(kr, model)
    has_scaling = JutulDarcy.endpoint_scaling_is_active(kr)
    # Check that scaling interface works
    @test has_scaling == (scaling != JutulDarcy.NoKrScale())
    if has_scaling
        for k in scalers
            @test k in deps
        end
        for k in setdiff(all_scalers, scalers)
            @test !(k in deps)
        end
        if should_use_hysteresis
            for k in scalers
                @test Symbol("$(k)i") in deps
            end
            for k in setdiff(all_scalers, scalers)
                @test !(Symbol("$(k)i") in deps)
            end
        end
    end
    # Check hyst
    has_hyst = JutulDarcy.hysteresis_is_active(kr)
    @test has_hyst == should_use_hysteresis
    if has_hyst
        @test :MaxSaturations in deps
    end
    has_connate = :ConnateWater in deps
    has_aqua = AqueousPhase() in JutulDarcy.get_phases(sys)
    if has_aqua
        if length(mph) == 3
            @test has_connate
        elseif has_scaling
            @test has_connate
        else
            @test !has_connate
        end
    else
        @test !has_connate
    end
    pv = pore_volume(domain)
    irate = 500*sum(pv)/tot_time
    bc = FlowBoundaryCondition(nc, p0/2)
    src1  = SourceTerm(1, irate, fractional_flow = f1)
    forces1 = setup_forces(model, sources = src1, bc = bc)

    src2  = SourceTerm(1, irate, fractional_flow = f2)
    forces2 = setup_forces(model, sources = src2, bc = bc)

    forces = typeof(forces1)[]
    for i in 1:nstep
        push!(forces, forces1)
    end
    for i in 1:nstep
        push!(forces, forces2)
    end
    # Make two different rates to get hysteresis effects
    parameters = setup_parameters(model, PhaseViscosities = mu)
    state0 = setup_state(model, Pressure = p0, Saturations = s0)
    states, report = simulate(state0, model, timesteps,
        forces = forces, parameters = parameters, info_level = -1)
    @test length(states) == length(timesteps)
end
##
@testset "ReservoirRelativePermeability" begin
    all_hyst = [
        JutulDarcy.NoHysteresis(),
        JutulDarcy.ImbibitionOnlyHysteresis(),
        JutulDarcy.KilloughHysteresis(),
        JutulDarcy.CarlsonHysteresis(),
        JutulDarcy.JargonHysteresis()
    ]
    all_endscale = [
        JutulDarcy.NoKrScale(),
        JutulDarcy.TwoPointKrScale(),
        JutulDarcy.ThreePointKrScale()
    ]
    all_combinations = [:wog, :wo, :wg, :og]

    for ph in all_combinations
        if startswith("$ph", 'w')
            wet_hyst = all_hyst
        else
            wet_hyst = [JutulDarcy.NoHysteresis()]
        end
        @testset "$ph" begin
            for non_wetting_hysteresis in all_hyst
                @testset "Non-Wetting=$non_wetting_hysteresis" begin
                    for wetting_hysteresis in wet_hyst
                        @testset "Wetting=$wetting_hysteresis" begin
                            for scaling in all_endscale
                                @testset "Scaling=$scaling" begin
                                    solve_1d_kr(ph,
                                        scaling = scaling,
                                        non_wetting_hysteresis = non_wetting_hysteresis,
                                        wetting_hysteresis = wetting_hysteresis
                                    )
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end;
