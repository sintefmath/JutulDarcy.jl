using CUDA, JutulDarcy, MultiComponentFlash, StaticArrays, Test

const FlashFD = JutulDarcy.ForwardDiff

function cubic_flash_derivatives_kernel!(out, eos, options, initial)
    if CUDA.threadIdx().x == 1
        D = FlashFD.Dual{Nothing, Float64, 4}
        p = D(1e6, FlashFD.Partials{4, Float64}((1e6, 0.0, 0.0, 0.0)))
        T = D(300.0, FlashFD.Partials{4, Float64}((0.0, 1.0, 0.0, 0.0)))
        z1 = D(0.6, FlashFD.Partials{4, Float64}((0.0, 0.0, 1.0, 0.0)))
        z2 = D(0.1, FlashFD.Partials{4, Float64}((0.0, 0.0, 0.0, 1.0)))
        z = SVector(z1, z2, one(D) - z1 - z2)
        result = JutulDarcy.immutable_flash_result(
            initial, options, eos, p, T, z, 0.0)
        ratio = result.vapor.mole_fractions[1]/result.liquid.mole_fractions[1]
        @inbounds for i in 1:4
            out[i] = result.V.partials[i]
            out[4+i] = ratio.partials[i]
        end
    end
    return nothing
end

function four_k_evaluator(cond)
    return SVector(
        0.05*(cond.p/1e6)*(cond.T/300)^0.2,
        5.0*(1e6/cond.p)*(cond.T/300)^(-0.1),
        0.2*(cond.T/300)^0.3,
        2.0*(cond.p/1e6)^0.2)
end

function kvalue_flash_derivatives_kernel!(out, eos, initial)
    if CUDA.threadIdx().x == 1
        D = FlashFD.Dual{Nothing, Float64, 5}
        p = D(1e6, FlashFD.Partials{5, Float64}((1e6, 0.0, 0.0, 0.0, 0.0)))
        T = D(300.0, FlashFD.Partials{5, Float64}((0.0, 1.0, 0.0, 0.0, 0.0)))
        z1 = D(0.2, FlashFD.Partials{5, Float64}((0.0, 0.0, 1.0, 0.0, 0.0)))
        z2 = D(0.3, FlashFD.Partials{5, Float64}((0.0, 0.0, 0.0, 1.0, 0.0)))
        z3 = D(0.2, FlashFD.Partials{5, Float64}((0.0, 0.0, 0.0, 0.0, 1.0)))
        z = SVector(z1, z2, z3, one(D) - z1 - z2 - z3)
        result = JutulDarcy.immutable_flash_result(
            initial, nothing, eos, p, T, z, 0.0)
        @inbounds for i in 1:5
            out[i] = result.V.partials[i]
        end
    end
    return nothing
end

@testset "Two-phase flash derivatives on CUDA" begin
    components = [
        MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224),
        MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011),
        MolecularProperty(0.0142, 2.10e6, 617.7, 6.098e-4, 0.488),
    ]
    cubic = make_eos_immutable(GenericCubicEOS(MultiComponentMixture(components)))
    method = SSIFlash()
    options = JutulDarcy.FlashResults{typeof(method), false, false}(
        method, false, 1e-8, 10.0, false, false)
    initial = JutulDarcy.static_flashed_mixture(
        cubic, FlashFD.Dual{Nothing, Float64, 4})
    output = CUDA.zeros(Float64, 8)
    CUDA.@sync CUDA.@cuda threads=1 cubic_flash_derivatives_kernel!(
        output, cubic, options, initial)
    derivatives = Array(output)
    @test all(isfinite, derivatives)
    @test derivatives[1:4] ≈ [-0.08319514018647238, 0.0013423289190622903,
        1.221240134768857, 1.3361248948491056] rtol=1e-5
    @test derivatives[5:8] ≈ [-4.245486438364429, 0.06981929786445129,
        -0.004528260740076506, 0.020547589014410228] rtol=1e-5

    mixture = MultiComponentMixture(
        ["CarbonDioxide", "Water", "Methane", "Ethane"])
    kvalue = make_eos_immutable(KValuesEOS(four_k_evaluator, mixture))
    initial_k = JutulDarcy.static_flashed_mixture(
        kvalue, FlashFD.Dual{Nothing, Float64, 5})
    koutput = CUDA.zeros(Float64, 5)
    CUDA.@sync CUDA.@cuda threads=1 kvalue_flash_derivatives_kernel!(
        koutput, kvalue, initial_k)
    kderivatives = Array(koutput)
    @test all(isfinite, kderivatives)
    @test kderivatives ≈ [-0.05201383595200247, 4.285767787159687e-5,
        -1.4715679450380024, 0.4180798543435838,
        -1.1980779158365582] rtol=1e-5
end
