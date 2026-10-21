using Test, Jutul, GeoEnergyIO, JutulDarcy

@testset "OPEN zero-rate deck controls" begin
    JD = JutulDarcy
    sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase(), VaporPhase());
        reference_densities = [1000.0, 800.0, 1.0])
    streams = JD.setup_well_streams()
    for kind in ("ORAT", "WRAT", "GRAT", "LRAT", "RESV"), is_hist in (false, true)
        ctrl, limits, flag = JD.producer_control(sys, "OPEN", kind,
            0.0, 0.0, 0.0, 0.0, 180e5; is_hist, resv = 0.0, factor = 0.8)
        @test ctrl isa ProducerControl
        @test ctrl.target isa TotalRateTarget
        @test iszero(ctrl.target.value)
        @test JD.is_zero_rate_control(ctrl)
        @test ctrl.factor == 0.8
        @test flag == "OPEN"
        @test haskey(limits, :bhp)
        if is_hist
            # Retain the history keyword's limit so WELTARG can update it
            # without leaving a zero total-rate cap after reopening.
            @test !haskey(limits, :rate)
            if kind != "RESV"
                @test haskey(limits, Symbol(lowercase(kind)))
            end
        end
    end
    for kind in ("RATE", "RESV"), is_hist in (false, true)
        ctrl, limits, flag = JD.injector_control(sys, streams, "INJ", "OPEN",
            "WATER", kind, 0.0, 0.0, 220e5; is_hist, factor = 0.7, temperature = 320.0)
        @test ctrl isa InjectorControl
        @test ctrl.target isa TotalRateTarget
        @test iszero(ctrl.target.value)
        @test JD.is_zero_rate_control(ctrl)
        @test ctrl.factor == 0.7
        @test ctrl.injection_mixture == [1.0, 0.0, 0.0]
        @test ctrl.mixture_density == 1000.0
        @test ctrl.temperature == 320.0
        @test flag == "OPEN"
        @test haskey(limits, :bhp)
    end
    for flag in ("SHUT", "STOP")
        ctrl, limits, = JD.producer_control(sys, flag, "ORAT", 0.0, 0.0, 0.0, 0.0, 180e5)
        @test ctrl isa DisabledControl
        @test isnothing(limits)
        ctrl, limits, = JD.injector_control(sys, streams, "INJ", flag,
            "WATER", "RATE", 0.0, 0.0, 220e5)
        @test ctrl isa DisabledControl
        @test isnothing(limits)
    end
    # Only the selected rate target determines closure. Zero secondary limits
    # must not close a BHP-controlled well or a well with a nonzero rate target.
    ctrl, = JD.producer_control(sys, "OPEN", "ORAT", 1.0, 0.0, 0.0, 0.0, 180e5)
    @test ctrl.target isa SurfaceOilRateTarget
    @test ctrl.target.value == -1.0
    ctrl, = JD.producer_control(sys, "OPEN", "BHP", 0.0, 0.0, 0.0, 0.0, 180e5)
    @test ctrl.target isa BottomHolePressureTarget
    ctrl, = JD.injector_control(sys, streams, "INJ", "OPEN", "WATER", "BHP", 0.0, 0.0, 220e5)
    @test ctrl.target isa BottomHolePressureTarget

    @testset "WELTARG zero and reopening" begin
        for initial in (ProducerControl(SurfaceOilRateTarget(-1.0)),
                InjectorControl(TotalRateTarget(1.0), [1.0, 0.0, 0.0]; density = 1000.0))
            keyword = initial isa ProducerControl ? "ORAT" : "RATE"
            controls = Dict{String, Any}("W" => initial)
            limits = Dict{String, Any}("W" => JD.default_limits(initial))
            JD.apply_weltarg!(controls, limits, ("W", keyword, 0.0))
            @test JD.is_zero_rate_control(controls["W"])
            JD.apply_weltarg!(controls, limits, ("W", "BHP", 200e5))
            @test JD.is_zero_rate_control(controls["W"])
            JD.apply_weltarg!(controls, limits, ("W", keyword, 2.0))
            @test typeof(controls["W"]) == typeof(initial)
            @test abs(controls["W"].target.value) == 2.0
        end
    end
end

@testset "SPE1 OPEN ORAT zero from DATA file" begin
    spe1 = read(GeoEnergyIO.test_input_file_path("SPE1", "SPE1.DATA"), String)
    pattern = r"'PROD'\s+'OPEN'\s+'ORAT'\s+20000"
    @test occursin(pattern, spe1)
    zero_rate_deck = replace(spe1, pattern => "'PROD' 'OPEN' 'ORAT' 0.0")
    mktempdir() do dir
        path = joinpath(dir, "Z.DATA")
        write(path, zero_rate_deck)
        case = setup_case_from_data_file(path)
        ctrl = case.forces[1][:Facility].control[:PROD]
        @test ctrl isa ProducerControl
        @test ctrl.target isa TotalRateTarget
        @test iszero(ctrl.target.value)
        @test JutulDarcy.is_zero_rate_control(ctrl)
        # OPEN zero rate retains active reservoir connections.
        mask = case.forces[1][:PROD].mask
        @test isnothing(mask) || any(>(0), mask.values)
        result = simulate_reservoir(case[1:2]; linear_solver = :fgmres, info_level = -1)
        @test last(result.time) ≈ sum(case.dt[1:2])
        @test all(iszero, result.wells[:PROD][:mass_rate])
        @test all(q -> abs(q) < 1e-9, result.wells[:PROD][:rate])
        @test all(isfinite, result.wells[:PROD][:bhp])
    end
end
