using Test, Jutul, JutulDarcy

function co2store_column_deck()
    return """
RUNSPEC
DIMENS
1 1 20 /
OIL
GAS
CO2STOR
METRIC
GRID
DXV
100 /
DYV
100 /
DZV
20*10 /
TOPS
798 /
PORO
20*0.3 /
PERMX
20*2000 /
COPY
PERMX PERMY /
PERMX PERMZ /
/
PROPS
ROCK
10 1.6e-6 /
SGOF
0.0 0.0 1.0 0.025
0.1 0.0 0.740 0.026
0.2 0.009 0.528 0.027
0.3 0.030 0.359 0.028
0.4 0.070 0.230 0.030
0.5 0.136 0.136 0.032
0.6 0.230 0.070 0.035
0.7 0.359 0.030 0.038
0.8 0.528 0.009 0.044
0.9 0.740 0.000 0.057 /
SALINITY
0.7 /
SOLUTION
EQUIL
863 146.5583 5050 0 100 0 1 0 0 /
/
RSVD
800 0
4150 0 /
RTEMPVD
800 32.7134
1000 41 /
SCHEDULE
TSTEP
1 /
"""
end

@testset "OPM CO2STORE conversion" begin
    mktempdir() do dir
        path = joinpath(dir, "COLUMN.DATA")
        write(path, co2store_column_deck())
        data = JutulDarcy.GeoEnergyIO.parse_data_file(path; silent = true)
        original = deepcopy(data)
        convert_deck(d) = JutulDarcy.convert_co2store_to_co2_brine(d; verbose = false)
        converted = convert_deck(data)
        rs = converted["RUNSPEC"]
        @test isequal(data, original)
        @test !haskey(rs, "CO2STOR")
        @test rs["JUTUL_CO2BRINE_PHYSICS"] == :immiscible
        @test only(converted["SOLUTION"]["EQUIL"])[3:4] == [100.0, 0.0]
        @test rs["SALTS"] == ["NaCl"]
        salt_fraction = 0.7/(0.7 + 1/0.01801528)
        @test rs["SALT_MOLE_FRACTIONS"] ≈ [salt_fraction]

        # Both keyword spellings and preconverted dictionaries use the same path.
        alias = deepcopy(data)
        delete!(alias["RUNSPEC"], "CO2STOR")
        alias["RUNSPEC"]["CO2STORE"] = true
        @test isequal(convert_deck(alias), converted)
        @test convert_deck(converted) === converted

        # Commercial compositional inputs retain their existing mixing/contact
        # convention rather than inheriting OPM's black-oil keyword defaults.
        compositional = deepcopy(data)
        compositional["RUNSPEC"]["COMPS"] = 2
        compositional = convert_deck(compositional)
        @test !haskey(compositional["RUNSPEC"], "JUTUL_CO2BRINE_PHYSICS")
        @test only(compositional["SOLUTION"]["EQUIL"])[3] == 5050.0

        for keyword in ("DISGAS", "DISGASW", "VAPOIL", "VAPWAT")
            mixing = deepcopy(data)
            mixing["RUNSPEC"][keyword] = true
            @test convert_deck(mixing)["RUNSPEC"]["JUTUL_CO2BRINE_PHYSICS"] == :kvalue
        end

        case = setup_case_from_data_file(data)
        rmodel = reservoir_model(case.model)
        state0 = case.state0[:Reservoir]
        @test rmodel.system isa ImmiscibleSystem
        @test all(isone, state0[:Saturations][1, :])
        @test all(iszero, state0[:Saturations][2, :])
        depths = rmodel.data_domain[:cell_centroids][3, :]
        datum_cell = argmin(abs.(depths .- 863.0))
        @test depths[datum_cell] ≈ 863.0
        @test state0[:Pressure][datum_cell] ≈ 146.5583e5 atol = 1.0
        # Salinity must reach the PVT model, not merely survive in the dictionary.
        density = rmodel.secondary_variables[:PhaseMassDensities].tab(146.5583e5, 310.0)[1]
        freshwater = JutulDarcy.CO2Properties.compute_co2_brine_props(146.5583e5, 310.0)[:density][1]
        @test density > freshwater + 20.0

        result = simulate_reservoir(case; info_level = -1, linear_solver = nothing,
            timesteps = :none, error_on_incomplete = true)
        @test length(result.states) == 1
        # A closed, initially hydrostatic column should not generate a transient.
        @test maximum(abs.(result.states[end][:Pressure] - state0[:Pressure])) < 100.0
        @test all(iszero, result.states[end][:Saturations][2, :])

        # Explicit user physics settings take precedence over conversion defaults.
        override = setup_case_from_data_file(converted; co2_physics = :kvalue)
        @test reservoir_model(override.model).system isa JutulDarcy.CompositionalSystem
        live = deepcopy(data)
        live["RUNSPEC"]["DISGAS"] = true
        live["RUNSPEC"]["VAPOIL"] = true
        @test reservoir_model(setup_case_from_data_file(live).model).system isa JutulDarcy.CompositionalSystem
    end
end
