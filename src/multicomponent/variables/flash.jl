@enum PERFORMED_FLASH_TYPE FLASH_SINGLE_PHASE_TESTED FLASH_SINGLE_PHASE_BYPASSED FLASH_FULL FLASH_RESTARTED FLASH_PURE_WATER

struct FlashResults{M, StabilityBypass, ReuseGuess} <: ScalarVariable
    method::M
    use_threads::Bool
    tolerance::Float64
    tolerance_bypass::Float64
    stability_bypass::Bool
    reuse_guess::Bool
end

function FlashResults(model::SimulationModel;
        method::MultiComponentFlash.AbstractFlash = SSIFlash(),
        threads = Threads.nthreads() > 1,
        tolerance = 1e-8,
        tolerance_bypass = 10,
        reuse_guess = false,
        stability_bypass = false,
        kwarg...)
    method isa SSIFlash || throw(ArgumentError(
        "The immutable compositional flash currently supports SSIFlash only."))
    return FlashResults{typeof(method), stability_bypass, reuse_guess}(
        method, threads, Float64(tolerance), Float64(tolerance_bypass),
        stability_bypass, reuse_guess)
end

@inline function static_flashed_mixture(
        eos::MultiComponentFlash.AbstractEOS, ::Type{T}) where T
    n = MultiComponentFlash.number_of_components(eos)
    K = SVector{n, Float64}(ntuple(_ -> NaN, n))
    xy = SVector{n, T}(ntuple(_ -> zero(T), n))
    cond = (p = NaN, T = NaN, z = K)
    return FlashedMixture2Phase(
        MultiComponentFlash.unknown_phase_state_lv,
        K, zero(T), xy, xy, zero(T), zero(T), NaN, cond)
end

function default_value(model, ::FlashResults,
        T = Jutul.float_type(model.context))
    return static_flashed_mixture(model.system.equation_of_state, T)
end

function initialize_variable_value(model, pvar::FlashResults,
        val::AbstractDict; need_value = false,
        T = Jutul.float_type(model.context))
    @assert need_value == false
    n = number_of_entities(model, pvar)
    return fill(default_value(model, pvar, T), n)
end

function Jutul.initialize_variable_value!(state, model, pvar::FlashResults,
        symb, val::AbstractDict; kwarg...)
    state[symb] = initialize_variable_value(model, pvar, val; kwarg...)
end

function initialize_variable_ad!(state, model, pvar::FlashResults, symb,
        npartials, diag_pos; context = DefaultContext(), kwarg...)
    n = number_of_entities(model, pvar)
    sample = get_ad_entity_scalar(1.0, npartials, diag_pos; kwarg...)
    state[symb] = fill(
        static_flashed_mixture(model.system.equation_of_state, typeof(sample)), n)
    return state
end

@jutul_secondary function update_flash!(flash_results, fr::FlashResults,
        model, Pressure, Temperature, OverallMoleFractions, ix)
    eos = model.system.equation_of_state
    flash_entity_loop!(flash_results, fr, eos, Pressure, Temperature,
        OverallMoleFractions, nothing, ix)
end

@jutul_secondary function update_flash!(flash_results, fr::FlashResults,
        model::LVCompositionalModel3Phase, Pressure, Temperature,
        OverallMoleFractions, ImmiscibleSaturation, ix)
    eos = model.system.equation_of_state
    flash_entity_loop!(flash_results, fr, eos, Pressure, Temperature,
        OverallMoleFractions, ImmiscibleSaturation, ix)
end

@inline immiscible_saturation(::Nothing, i) = 0.0
@inline immiscible_saturation(saturation, i) = @inbounds saturation[i]

@inline function flash_entity_loop!(flash_results, fr, eos, Pressure,
        Temperature, OverallMoleFractions, sw, ix)
    @inbounds for i in ix
        old = flash_results[i]
        flash_results[i] = immutable_flash_result(old, fr, eos,
            Pressure[i], Temperature[i], OverallMoleFractions,
            immiscible_saturation(sw, i), i)
    end
    return flash_results
end

@inline function cell_composition(::Val{N}, composition, cell) where N
    T = eltype(composition)
    return SVector{N, T}(ntuple(i -> @inbounds(composition[i, cell]), Val(N)))
end

@inline function cell_composition(::Val{N}, composition::AbstractVector, cell) where N
    T = eltype(composition)
    return SVector{N, T}(ntuple(i -> @inbounds(composition[i]), Val(N)))
end

@inline function numeric_composition(z::SVector{N}) where N
    return SVector{N, Float64}(ntuple(Val(N)) do i
        max(Float64(compositional_primal(z[i])),
            MultiComponentFlash.MINIMUM_COMPOSITION)
    end)
end

@inline compositional_primal(x) = x
@inline compositional_primal(x::ForwardDiff.Dual) = ForwardDiff.value(x)

@inline previous_stability_storage(f, ::Val{false}) = nothing
@inline function previous_stability_storage(f, ::Val{true})
    return MultiComponentFlash.StaticStabilityStorage(
        f.flash_cond, f.critical_distance)
end

@inline function full_numeric_flash(f,
        fr::FlashResults{M, StabilityBypass}, eos, cond) where {
        M, StabilityBypass}
    config = MultiComponentFlash.StaticConfig()
    K0 = initial_guess_K(eos, cond, config)
    storage = previous_stability_storage(f, Val(StabilityBypass))
    V, K, report = flash_2ph!(config, K0, eos, cond, NaN;
        method = SSIFlash(),
        extra_out = true,
        tolerance = fr.tolerance,
        z_min = nothing,
        stability_storage = storage,
        stability_bypass = StabilityBypass,
        bypass_tolerance = fr.tolerance_bypass,
        check = false,
        verbose = false)
    return V, K, report.stability_result
end

@inline numeric_flash(f, fr::FlashResults{M, StabilityBypass, false},
    eos, cond) where {M, StabilityBypass} =
    full_numeric_flash(f, fr, eos, cond)

@inline function numeric_flash(f,
        fr::FlashResults{M, StabilityBypass, true}, eos, cond) where {
        M, StabilityBypass}
    if f.state == MultiComponentFlash.two_phase_lv
        config = MultiComponentFlash.StaticConfig()
        V0 = Float64(compositional_primal(f.V))
        V, K, report = flash_2ph!(config, f.K, eos, cond, V0;
            method = SSIFlash(),
            maxiter = 20,
            extra_out = true,
            tolerance = fr.tolerance,
            z_min = nothing,
            stability_bypass = false,
            check = false,
            verbose = false)
        valid = report.converged && isfinite(V) &&
            1e-6 < V < 1.0 - 1e-6
        trivial = true
        @inbounds for K_i in K
            trivial &= abs(K_i - 1.0) < 1e-6
        end
        if valid && !trivial
            return V, K, report.stability_result
        end
    end
    return full_numeric_flash(f, fr, eos, cond)
end

@inline function pure_immiscible_flash(f, eos, cond)
    config = MultiComponentFlash.StaticConfig()
    K = initial_guess_K(eos, cond, config)
    report = MultiComponentFlash.StabilityReport(
        stable_liquid = true, stable_vapor = true)
    storage = MultiComponentFlash.StaticStabilityStorage(cond, NaN)
    stability = MultiComponentFlash.StaticStabilityResult(
        true, report, K, storage, false)
    return NaN, K, stability
end

@inline function equilibrium_ad(eos, cond, vapor_fraction, K_numeric,
        pressure::AbstractFloat, tolerance)
    T = typeof(cond.p)
    return convert(T, vapor_fraction), SVector{length(K_numeric), T}(K_numeric)
end

@inline function equilibrium_ad(eos, cond, vapor_fraction, K_numeric,
        pressure::ForwardDiff.Dual, tolerance)
    T = typeof(cond.p)
    config = MultiComponentFlash.StaticConfig()
    K0 = initial_guess_K(eos, cond, config)
    V, K, _ = flash_2ph!(config, K0, eos, cond,
        convert(T, vapor_fraction);
        method = SSIFlash(),
        extra_out = true,
        tolerance = tolerance,
        z_min = nothing,
        stability_bypass = false,
        check = false,
        verbose = false)
    return V, K
end

@inline function phase_compressibility(eos, cond, phase)
    T = typeof(cond.p)
    phase_cond = (p = cond.p, T = cond.T, z = cond.z, phase = phase)
    forces = MultiComponentFlash.static_force_coefficients(
        eos, phase_cond, T)
    scalars = force_scalars(eos, phase_cond, forces)
    return mixture_compressibility_factor(
        eos, phase_cond, forces, scalars)
end


@generated function phase_mole_fractions(z::SVector{N, T},
        K::SVector{N, T}, V::T) where {N, T}
    x_values = [:(liquid_mole_fraction(z[$i], K[$i], V)) for i in 1:N]
    y_values = [:(vapor_mole_fraction(x[$i], K[$i])) for i in 1:N]
    return quote
        x = SVector{N, T}(($(x_values...),))
        y = SVector{N, T}(($(y_values...),))
        (x, y)
    end
end

@inline function single_phase_flash_result(eos, cond, cond_numeric,
        K_numeric, stability, is_vapor::Bool)
    Num = typeof(cond.p)
    state = is_vapor ? MultiComponentFlash.single_phase_v :
        MultiComponentFlash.single_phase_l
    V = convert(Num, is_vapor)
    x = y = cond.z
    if is_vapor
        Z_l = Z_v = phase_compressibility(eos, cond, Val(:vapor))
    else
        Z_l = Z_v = phase_compressibility(eos, cond, Val(:liquid))
    end
    return state, V, x, y, Z_l, Z_v
end

@inline function two_phase_flash_result(eos, cond, P, V_numeric,
        K_numeric, tolerance)
    V, K = equilibrium_ad(
        eos, cond, V_numeric, K_numeric, P, tolerance)
    x, y = phase_mole_fractions(cond.z, K, V)
    liquid = (p = cond.p, T = cond.T, z = x)
    vapor = (p = cond.p, T = cond.T, z = y)
    Z_l = phase_compressibility(eos, liquid, Val(:liquid))
    Z_v = phase_compressibility(eos, vapor, Val(:vapor))
    return MultiComponentFlash.two_phase_lv, V, x, y, Z_l, Z_v
end

@inline function immutable_flash_result(f,
        fr::FlashResults{M, StabilityBypass},
        eos::GenericCubicEOS{E, R, N}, P, temperature,
        OverallMoleFractions, Sw, cell = 1) where {
        M, StabilityBypass, E, R, N}
    z = cell_composition(Val(N), OverallMoleFractions, cell)
    z_numeric = numeric_composition(z)
    cond_numeric = (
        p = Float64(compositional_primal(P)),
        T = Float64(compositional_primal(temperature)),
        z = z_numeric)
    if is_pure_single_phase(compositional_primal(Sw))
        V_numeric, K_numeric, stability =
            pure_immiscible_flash(f, eos, cond_numeric)
    else
        V_numeric, K_numeric, stability =
            numeric_flash(f, fr, eos, cond_numeric)
    end

    Num = typeof(P + temperature + first(z))
    cond = (p = convert(Num, P), T = convert(Num, temperature), z = z)
    if isnan(V_numeric)
        is_vapor = single_phase_label(eos, cond_numeric) > 0.5
        state, V, x, y, Z_l, Z_v = single_phase_flash_result(
            eos, cond, cond_numeric, K_numeric, stability, is_vapor)
    else
        state, V, x, y, Z_l, Z_v = two_phase_flash_result(
            eos, cond, P, V_numeric, K_numeric, fr.tolerance)
    end

    K_out = K_numeric
    if StabilityBypass
        critical_distance = Float64(
            compositional_primal(stability.storage.critical_distance))
        flash_cond = stability.storage.reference
    else
        critical_distance = NaN
        flash_cond = cond_numeric
    end
    return FlashedMixture2Phase(state, K_out, V, x, y, Z_l, Z_v,
        critical_distance, flash_cond, stability.report)
end

@inline function immutable_flash_result(f, fr,
        eos::KValuesEOS{E, R, N}, P, temperature,
        OverallMoleFractions, Sw, cell = 1) where {E, R, N}
    z = cell_composition(Val(N), OverallMoleFractions, cell)
    Num = Base.promote_type(typeof(P), typeof(temperature), eltype(z))
    z = SVector{N, Num}(z)
    cond = (p = convert(Num, P), T = convert(Num, temperature), z = z)
    K = SVector{N, Num}(initial_guess_K(eos, cond))
    V = MultiComponentFlash.solve_rachford_rice(K, z)
    if V <= zero(V)
        state = MultiComponentFlash.single_phase_l
        V = zero(V)
        x = y = z
    elseif V >= one(V)
        state = MultiComponentFlash.single_phase_v
        V = one(V)
        x = y = z
    else
        state = MultiComponentFlash.two_phase_lv
        x = SVector{N, Num}(ntuple(Val(N)) do component
            liquid_mole_fraction(z[component], K[component], V)
        end)
        y = SVector{N, Num}(ntuple(Val(N)) do component
            vapor_mole_fraction(x[component], K[component])
        end)
    end
    K_out = SVector{N, Float64}(ntuple(
        i -> Float64(value(K[i])), Val(N)))
    cond_numeric = (p = Float64(value(P)), T = Float64(value(temperature)),
        z = numeric_composition(z))
    return FlashedMixture2Phase(state, K_out, V, x, y,
        one(Num), one(Num), NaN, cond_numeric, f.flash_stability)
end

@inline replace_flash_value(old, new) = old - value(old) + value(new)

@inline function replace_flash_vector(old::SVector{N}, new) where N
    return SVector{N}(ntuple(
        i -> replace_flash_value(old[i], new[i]), Val(N)))
end

@inline function replace_flash_values(old::FlashedMixture2Phase, next)
    liquid_xy = replace_flash_vector(
        old.liquid.mole_fractions, next.liquid.mole_fractions)
    vapor_xy = replace_flash_vector(
        old.vapor.mole_fractions, next.vapor.mole_fractions)
    liquid = FlashedPhase(liquid_xy,
        replace_flash_value(old.liquid.Z, next.liquid.Z))
    vapor = FlashedPhase(vapor_xy,
        replace_flash_value(old.vapor.Z, next.vapor.Z))
    V = replace_flash_value(old.V, next.V)
    return FlashedMixture2Phase(next.state, next.K, V, liquid, vapor;
        vec_type = typeof(liquid_xy),
        critical_distance = next.critical_distance,
        cond = next.flash_cond,
        stability_report = next.flash_stability)
end

function Jutul.update_values!(vals::AbstractVector{<:FlashedMixture2Phase},
        next::AbstractVector{<:FlashedMixture2Phase})
    @inbounds for i in eachindex(vals, next)
        vals[i] = replace_flash_values(vals[i], next[i])
    end
    return vals
end


function Jutul.update_values!(vals::AbstractVector{<:FlashedMixture2Phase},
        next::AbstractVector{<:FlashedMixture2Phase},
        context::Jutul.KernelAbstractionsContext)
    next_backend = Adapt.adapt(context, next)
    function update(i)
        @inbounds vals[i] = replace_flash_values(vals[i], next_backend[i])
    end
    Jutul.threaded_loop_minbatch(update, length(vals), context)
    return vals
end
