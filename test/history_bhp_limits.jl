using Test, Jutul, GeoEnergyIO, JutulDarcy

@testset "History producer BHP limits" begin
    sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase(), VaporPhase());
        reference_densities = [1000.0, 800.0, 1.0])
    function parse_history_schedule(records)
        mktemp() do path, io
            write(io, """
RUNSPEC
METRIC
OIL
WATER
GAS
SCHEDULE
WELSPECS
 P G 1 1 1000 OIL /
/
$records
END
""")
            close(io)
            data = parse_data_file(path; silent = true)
            return JutulDarcy.parse_control_steps(data["RUNSPEC"], Dict(), data["SCHEDULE"], sys)
        end
    end
    history(mode; flag = "OPEN", bhp = 200) = "WCONHIST\n P $flag $mode 100 10 0 2* 0 $bhp /\n/\n"
    step = "TSTEP\n 1 /\n"
    bhp_limit(p) = [lim["P"].bhp/si_unit(:bar) for lim in p.limits]
    for mode in ("ORAT", "RESV")
        records = history(mode)*"WELTARG\n P BHP 50 /\n/\n"*step*
            history(mode; bhp = 250)*step
        parsed = parse_history_schedule(records)
        @test bhp_limit(parsed) ≈ [50, 50]
        @test parsed.controls[2]["P"] isa ProducerControl
        @test !(parsed.controls[2]["P"].target isa BottomHolePressureTarget)

        # Changing the field default must not erase an explicit history limit.
        parsed = parse_history_schedule(records*"FBHPDEF\n 75 /\n"*history(mode)*step)
        @test bhp_limit(parsed) ≈ [50, 50, 50]

        # SHUT records retain the history property for WELOPEN and later updates.
        parsed = parse_history_schedule(records*history(mode; flag = "SHUT")*step*
            "WELOPEN\n P OPEN /\n/\n"*step*history(mode)*step)
        @test bhp_limit(parsed) ≈ fill(50, 5)
        @test parsed.controls[3]["P"] isa DisabledControl
        @test parsed.controls[4]["P"] isa ProducerControl
    end

    # Observed BHP is not a rate-history pressure limit.
    parsed = parse_history_schedule(history("ORAT")*step)
    @test only(bhp_limit(parsed)) ≈ 1.01*si_unit(:atm)/si_unit(:bar)
    parsed = parse_history_schedule("FBHPDEF\n 75 /\n"*history("ORAT")*step)
    @test only(bhp_limit(parsed)) ≈ 75

    # Re-entering rate history from BHP, prediction, or injection resets the
    # limit to FBHPDEF when no explicit producer history limit was set.
    for previous in (history("BHP"),
            "WCONPROD\n P OPEN BHP 5* 200 /\n/\n",
            "WCONINJE\n P WATER OPEN RATE 100 1* 300 /\n/\n",
            "WCONINJH\n P WATER OPEN 100 300 /\n/\n")
        parsed = parse_history_schedule("FBHPDEF\n 75 /\n"*previous*step*history("ORAT")*step)
        @test last(bhp_limit(parsed)) ≈ 75
    end
    for injection in ("WCONINJE\n P WATER OPEN RATE 100 1* 300 /\n/\n",
            "WCONINJH\n P WATER OPEN 100 300 /\n/\n")
        parsed = parse_history_schedule(history("ORAT")*"WELTARG\n P BHP 50 /\n/\n"*step*
            injection*step*"FBHPDEF\n 75 /\n"*history("ORAT")*step)
        @test last(bhp_limit(parsed)) ≈ 50
    end
end
