using JutulDarcy, Test, LinearAlgebra

function test_mrst_case(casename; tol = 0.01, wtol = tol, otol = tol, gtol = tol, bhptol = tol)
    data_path = joinpath(@__DIR__, "mrst", "$(casename).mat")
    if isfile(data_path)
        result = simulate_mrst_case(data_path, info_level = -1, verbose = false)
        # ws = result.wells
        case = result.extra[:case]
        ref = result.extra[:mrst]["extra"][1]["mrst_solution"]
        dt = case.dt
        ws = full_well_outputs(case.model, result.result.states, case.forces)
        output = Dict();
        for k in keys(ws)
            output[k] = Dict()
            @testset "$k" begin
                ix = findfirst(isequal("$k"), vec(ref["names"]))
                w = (abs.(ws[k][:mass]) .> 1e-10).*dt
                for wfield in [:bhp, :wrat, :orat, :grat]
                    if haskey(ws[k], wfield)
                        is_bhp = false
                        if wfield == :bhp
                            mrst = ref["bhp"]
                            is_bhp = true
                            ϵ = bhptol
                        elseif wfield == :wrat
                            mrst = ref["qWs"]
                            ϵ = wtol
                        elseif wfield == :orat
                            mrst = ref["qOs"]
                            ϵ = otol
                        else
                            @assert wfield == :grat
                            mrst = ref["qGs"]
                            ϵ = gtol
                        end
                        jutul = ws[k][wfield]
                        m = mrst[:, ix]
                        output[k][wfield] = (jutul = jutul, mrst = m)
                        ref_norm = norm(m.*w, 1)
                        # Don't test pure noise. MRST and Jutul do not use exactly the same wells.
                        if ref_norm > 1e-8*sum(dt)
                            @testset "$wfield" begin
                                err = norm(jutul.*w - m.*w, 1)/ref_norm
                                if err > ϵ
                                    @info "$casename $k: $wfield" err ref_norm #jutul m ws[k][:mass] w
                                    println("Compare")
                                    display(hcat(jutul, m, w))
                                end
                                @test err < ϵ
                            end
                        end
                    end
                end
            end
        end
    else
        @warn "Did not find test case." data_path
    end
    return output
end
##
@testset "SPE1" begin
    test_mrst_case("spe1")
end
##
@testset "SPE3" begin
    test_mrst_case("spe3", wtol = 0.2, gtol = 0.02)
end
##
@testset "SPE9" begin
    test_mrst_case("spe9", wtol = Inf, gtol = 0.1, otol = 0.05)
end
##
@testset "Egg" begin
    test_mrst_case("egg", otol = 0.02)
end
