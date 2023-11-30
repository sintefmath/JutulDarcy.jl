using Jutul, JutulDarcy, Test

function scalarization_bl_test_sim(;nc = 10, nph = 2)
    domain = get_1d_reservoir(nc)
    if nph == 2
        phases = (LiquidPhase(), VaporPhase())
    else
        @assert nph == 3
        phases = (LiquidPhase(), AqueousPhase(), VaporPhase())
    end
    sys = ImmiscibleSystem(phases)
    model = SimulationModel(domain, sys)
    sat = zeros(nph, nc)
    for i in 1:nc
        v = rand()
        sat[1, i] = v
        for j in 2:nph
            sat[j, i] = (1.0 - v)/(nph-1)
        end
    end
    state0 = setup_state(model, Pressure = rand(nc).*1e5, Saturations = sat)
    sim = Simulator(model, state0 = state0)
    return sim
end

@testset "scalarization 2ph" begin
    for case in (2, 3, "spe1")
        if case isa Int
            test_name = "$case phases"
            sim = scalarization_bl_test_sim(nph = case);
        else
            test_name = "spe1"
            tmp, = setup_case_from_mrst(case)
            sim = Simulator(tmp.model[:Reservoir], state0 = tmp.state0[:Reservoir])
        end
        @testset "$test_name" begin
            state = sim.storage.primary_variables
            model = sim.model
            pvar_vec = Jutul.scalarize_primary_variables(model, state)
            state_new = Dict()
            for (k, v) in pairs(state)
                new_v = copy(v)
                for i in eachindex(new_v)
                    new_v[i] = Jutul.replace_value(new_v[i], zero(new_v[i]))
                end
                state_new[k] = new_v
            end
            state_new = (; pairs(state_new)...)
            state_new = Jutul.descalarize_primary_variables!(state_new, model, pvar_vec)

            for (k_var, var_new) in pairs(state_new)
                var_old = state[k_var]
                @test size(var_old) == size(var_new)
                @test typeof(var_old) == typeof(var_new)
                for i in eachindex(var_old)
                    @test var_old[i] == var_new[i]
                    if eltype(var_old)<:JutulDarcy.BlackOilX
                        @test var_old[i].val.partials == var_new[i].val.partials
                    else
                        @test var_old[i].partials == var_new[i].partials
                    end
                end
            end
        end
    end
end
