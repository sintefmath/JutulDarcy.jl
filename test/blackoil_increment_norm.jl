using Jutul, JutulDarcy, Test, Random, LinearAlgebra
using Jutul: ForwardDiff
const JDIncrementTests = JutulDarcy

@testset "Fused Float64 black-oil increment diagnostics" begin
    Random.seed!(1234)
    jd = JDIncrementTests
    phases = (jd.OilAndGas, jd.OilOnly, jd.GasOnly)
    ctx = KernelAbstractionsContext(Jutul.KernelExecution.KernelAbstractions.CPU())
    model = (context=ctx,)
    dual = ForwardDiff.Dual{Nothing}(1.0,0.0,0.0,0.0)
    for n in (0,1,3,67,16385), stride in (1,3), phase in (nothing,phases...), val in (1.0,dual)
        storage = randn(stride,n)
        increments = view(storage,stride,:)
        X = [jd.BlackOilX(val,isnothing(phase) ? phases[mod1(i,3)] : phase,false) for i in 1:n]
        for last_value in (0.0,NaN,Inf)
            n > 0 && (increments[end] = last_value)
            reference = Jutul.increment_norm(increments,missing,missing,X,jd.BlackOilUnknown())
            fused = jd.cpu_blackoil_increment_norms(increments,X)
            for dX in (increments,transpose(reshape(increments,n,1)),adjoint(reshape(increments,n,1)))
                actual = Jutul.increment_norm(dX,missing,model,X,jd.BlackOilUnknown())
                for key in (:sg_sum,:rs_sum,:rv_sum)
                    @test isapprox(reference[key],fused[key];rtol=5e-14,atol=1e-14) || isequal(reference[key],fused[key])
                    @test isapprox(reference[key],actual[key];rtol=5e-14,atol=1e-14) || isequal(reference[key],actual[key])
                end
                for key in (:sg_max,:rs_max,:rv_max)
                    @test isequal(reference[key],fused[key])
                    @test isequal(reference[key],actual[key])
                end
            end
        end
    end
    n = 20000
    X = [jd.BlackOilX(dual,phases[mod1(i,3)],false) for i in 1:n]
    increments = randn(n)
    Jutul.increment_norm(increments,missing,model,X,jd.BlackOilUnknown())
    @test (@allocated Jutul.increment_norm(increments,missing,model,X,jd.BlackOilUnknown())) < 4096
    @test_throws DimensionMismatch jd.cpu_blackoil_increment_norms(increments,X[1:end-1])
    # Cell-major primary-variable views select a one-row matrix with a range,
    # then vec wraps that noncontiguous view in a ReshapedArray.
    buffer = randn(3n+6)
    block_view = reshape(view(buffer,4:3n+3),3,n)
    cell_major_dx = view(block_view,3:3,:)
    reference = Jutul.increment_norm(vec(cell_major_dx),missing,missing,X,jd.BlackOilUnknown())
    actual = Jutul.increment_norm(cell_major_dx,missing,model,X,jd.BlackOilUnknown())
    @test all(isapprox(reference[k],actual[k];rtol=5e-14) for k in (:sg_sum,:rs_sum,:rv_sum))
    @test all(isequal(reference[k],actual[k]) for k in (:sg_max,:rs_max,:rv_max))
    @test (@allocated Jutul.increment_norm(cell_major_dx,missing,model,X,jd.BlackOilUnknown())) < 4096
    # Narrower storage and limited contexts retain the existing reductions.
    dx32 = Float32.(increments)
    expected32 = Jutul.increment_norm(dx32,missing,missing,X,jd.BlackOilUnknown())
    @test Jutul.increment_norm(dx32,missing,model,X,jd.BlackOilUnknown()) == expected32
    if Threads.nthreads(:default) > 1
        limited = (context=KernelAbstractionsContext(Jutul.KernelExecution.KernelAbstractions.CPU();nthreads=1),)
        reference = Jutul.increment_norm(increments,missing,missing,X,jd.BlackOilUnknown())
        @test Jutul.increment_norm(increments,missing,limited,X,jd.BlackOilUnknown()) == reference
        @test (@allocated Jutul.increment_norm(increments,missing,limited,X,jd.BlackOilUnknown())) > 8n
    end
end
