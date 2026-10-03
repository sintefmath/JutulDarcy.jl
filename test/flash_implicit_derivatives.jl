using Jutul, JutulDarcy, MultiComponentFlash, StaticArrays, Test

const FD = JutulDarcy.ForwardDiff
const SCT = Jutul.SCT

function finite_difference(f, u, i, h)
    plus = copy(u)
    minus = copy(u)
    plus[i] += h
    minus[i] -= h
    return (f(plus) - f(minus))/(2h)
end

function check_two_phase_flash_derivatives(flash, u, steps;
        atol = 1e-7, rtol = 2e-4)
    result = flash(u)
    @test result.state == MultiComponentFlash.two_phase_lv
    outputs = (
        r -> r.V,
        r -> r.liquid.mole_fractions[1],
        r -> r.vapor.mole_fractions[1],
        r -> r.vapor.mole_fractions[1]/r.liquid.mole_fractions[1],
        r -> r.liquid.Z,
        r -> r.vapor.Z,
    )
    for output in outputs
        f = v -> output(flash(v))
        # A combined seed exercises the normal reservoir Jacobian path.
        gradient = FD.gradient(f, u)
        for i in eachindex(u)
            expected = finite_difference(f, u, i, steps[i])
            @test isapprox(gradient[i], expected; atol, rtol)
            # Single-input seeds catch a missing T or z derivative when p
            # itself is numeric.
            derivative = FD.derivative(t -> begin
                v = [j == i ? t : u[j] for j in eachindex(u)]
                f(v)
            end, u[i])
            @test isapprox(derivative, expected; atol, rtol)
        end
    end
end

@testset "Cubic flash implicit derivatives" begin
    components = [
        MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224),
        MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011),
        MolecularProperty(0.0142, 2.10e6, 617.7, 6.098e-4, 0.488),
    ]
    eos = GenericCubicEOS(MultiComponentMixture(components))
    method = SSIFlash()
    flash_options = JutulDarcy.FlashResults{typeof(method), false, false}(
        method, false, 1e-8, 10.0, false, false)
    initial = JutulDarcy.static_flashed_mixture(eos, Float64)
    flash = u -> JutulDarcy.immutable_flash_result(
        initial, flash_options, eos, 1e6*u[1], u[2],
        SVector(u[3], u[4], 1-u[3]-u[4]), 0.0)
    check_two_phase_flash_derivatives(
        flash, [1.0, 300.0, 0.6, 0.1], [1e-4, 1e-2, 1e-5, 1e-5])

    # Exercise the helper with heterogeneous condition types. The pressure
    # stays numeric while either temperature or composition carries AD.
    static_eos = make_eos_immutable(eos)
    numeric_cond = (p = 1e6, T = 300.0, z = @SVector [0.6, 0.1, 0.3])
    V_numeric, K_numeric = flash_2ph_immutable(static_eos, numeric_cond;
        check = false)
    mixed_result(cond) = JutulDarcy.two_phase_flash_result(
        static_eos, numeric_cond, cond, V_numeric, K_numeric)[2]
    temperature_derivative = FD.derivative(300.0) do T
        mixed_result((p = numeric_cond.p, T = T, z = numeric_cond.z))
    end
    composition_derivative = FD.derivative(0.6) do z1
        mixed_result((p = numeric_cond.p, T = numeric_cond.T,
            z = SVector(z1, 0.1, 0.9-z1)))
    end
    @test temperature_derivative ≈ FD.derivative(
        T -> flash([1.0, T, 0.6, 0.1]).V, 300.0) rtol=1e-5
    @test composition_derivative ≈ FD.derivative(
        z1 -> flash([1.0, 300.0, z1, 0.1]).V, 0.6) rtol=1e-5

    tracer_type = SCT.Dual{Float64,
        SCT.gradient_tracer_type(Set{Int})}
    traced = SCT.trace_input(tracer_type, [1.0, 300.0, 0.6, 0.1])
    result = flash(traced)
    @test Set(SCT.gradient(result.V)) == Set(1:4)
    @test Set(SCT.gradient(result.vapor.mole_fractions[1] /
        result.liquid.mole_fractions[1])) == Set(1:4)
end

@testset "K-value flash derivatives" begin
    mixture2 = MultiComponentMixture(["CarbonDioxide", "Water"])
    eos2 = KValuesEOS(cond -> SVector(
        0.05*(cond.p/1e6)*(cond.T/300)^0.2,
        5.0*(1e6/cond.p)*(cond.T/300)^(-0.1)), mixture2)
    initial2 = JutulDarcy.static_flashed_mixture(eos2, Float64)
    flash2 = u -> JutulDarcy.immutable_flash_result(
        initial2, nothing, eos2, 1e6*u[1], u[2],
        SVector(u[3], 1-u[3]), 0.0)
    check_two_phase_flash_derivatives(
        flash2, [1.0, 300.0, 0.4], [1e-4, 1e-2, 1e-5])

    mixture4 = MultiComponentMixture(
        ["CarbonDioxide", "Water", "Methane", "Ethane"])
    eos4 = KValuesEOS(cond -> SVector(
        0.05*(cond.p/1e6)*(cond.T/300)^0.2,
        5.0*(1e6/cond.p)*(cond.T/300)^(-0.1),
        0.2*(cond.T/300)^0.3,
        2.0*(cond.p/1e6)^0.2), mixture4)
    initial4 = JutulDarcy.static_flashed_mixture(eos4, Float64)
    flash4 = u -> JutulDarcy.immutable_flash_result(
        initial4, nothing, eos4, 1e6*u[1], u[2],
        SVector(u[3], u[4], u[5], 1-u[3]-u[4]-u[5]), 0.0)
    check_two_phase_flash_derivatives(
        flash4, [1.0, 300.0, 0.2, 0.3, 0.2],
        [1e-4, 1e-2, 1e-5, 1e-5, 1e-5])
end

@testset "Float32 cubic and K-value flash derivatives" begin
    components = [
        MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224),
        MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011),
        MolecularProperty(0.0142, 2.10e6, 617.7, 6.098e-4, 0.488),
    ]
    cubic = GenericCubicEOS(MultiComponentMixture(components))
    method = SSIFlash()
    options = JutulDarcy.FlashResults{typeof(method), false, false}(
        method, false, 1e-8, 10.0, false, false)
    initial_cubic = JutulDarcy.static_flashed_mixture(cubic, Float32)
    cubic_flash = u -> JutulDarcy.immutable_flash_result(
        initial_cubic, options, cubic, 1f6*u[1], u[2],
        SVector(u[3], u[4], 1f0-u[3]-u[4]), 0f0)

    mixture = MultiComponentMixture(["CarbonDioxide", "Water"])
    kvalue = KValuesEOS(cond -> SVector(
        0.05f0*(cond.p/1f6), 5f0*(1f6/cond.p)), mixture)
    initial_kvalue = JutulDarcy.static_flashed_mixture(kvalue, Float32)
    kvalue_flash = u -> JutulDarcy.immutable_flash_result(
        initial_kvalue, nothing, kvalue, 1f6*u[1], u[2],
        SVector(u[3], 1f0-u[3]), 0f0)

    for (flash, u, i) in ((cubic_flash, Float32[1, 300, 0.6, 0.1], 2),
            (kvalue_flash, Float32[1, 300, 0.4], 3))
        result = flash(u)
        @test result.state == MultiComponentFlash.two_phase_lv
        derivative = FD.derivative(t -> begin
            v = [j == i ? t : u[j] for j in eachindex(u)]
            flash(v).V
        end, u[i])
        @test derivative isa Float32
        @test isfinite(derivative)
        @test abs(derivative) > 1f-5
    end
end
