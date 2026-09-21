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
    # Specialize the two independent switches for allocation-free GPU kernels.
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

@inline function cell_composition(::Val{N}, composition, cell) where N
    T = eltype(composition)
    return SVector{N, T}(ntuple(i -> @inbounds(composition[i, cell]), Val(N)))
end

@inline function cell_composition(::Val{N}, composition::AbstractVector, cell) where N
    T = eltype(composition)
    return SVector{N, T}(ntuple(i -> @inbounds(composition[i]), Val(N)))
end

@inline function numeric_composition(z::SVector{N}, R = Float64) where N
    return SVector{N, R}(ntuple(Val(N)) do i
        max(R(compositional_primal(z[i])),
            MultiComponentFlash.MINIMUM_COMPOSITION)
    end)
end

@inline compositional_primal(x) = x
@inline compositional_primal(x::ForwardDiff.Dual) = ForwardDiff.value(x)

@inline function numeric_values(v::SVector{N}, R = Float64) where N
    return SVector{N, Float64}(ntuple(Val(N)) do i
        R(compositional_primal(v[i]))
    end)
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


@generated function phase_mole_fractions(z::SVector{N, Tz},
        K::SVector{N, Tk}, V::Tv) where {N, Tz, Tk, Tv}
    T = promote_type(Tz, Tk, Tv)
    x_values = [:(liquid_mole_fraction(z[$i], K[$i], V)) for i in 1:N]
    y_values = [:(vapor_mole_fraction(x[$i], K[$i])) for i in 1:N]
    return quote
        x = SVector{N, $T}(($(x_values...),))
        y = SVector{N, $T}(($(y_values...),))
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
