mutable struct InPlaceFlashBuffer
    z::Vector{Float64}
    forces
    function InPlaceFlashBuffer(n)
        z = zeros(n)
        new(z, nothing)
    end
end

@enum PERFORMED_FLASH_TYPE FLASH_SINGLE_PHASE_TESTED FLASH_SINGLE_PHASE_BYPASSED FLASH_FULL FLASH_RESTARTED FLASH_PURE_WATER
struct FlashResults{F_t, S_t, B_t} <: ScalarVariable
    storage::S_t
    method::F_t
    update_buffer::B_t
    use_threads::Bool
    tolerance::Float64
    tolerance_bypass::Float64
    stability_bypass::Bool
    reuse_guess::Bool
    function FlashResults(model::SimulationModel;
            method::MultiComponentFlash.AbstractFlash = SSIFlash(),
            threads = Threads.nthreads() > 1,
            n = degrees_of_freedom_per_entity(model, Cells()),
            tolerance = 1e-8,
            tolerance_bypass = 10,
            reuse_guess = false,
            stability_bypass = false
        )
        system = flow_system(model.system)
        eos = system.equation_of_state
        nc = MultiComponentFlash.number_of_components(eos)
        storage = []
        buffers = []
        if threads
            N = Threads.nthreads()
        else
            N = 1
        end
        for i = 1:N
            s = flash_storage(eos,
                method = method,
                inc_jac = true,
                inc_bypass = stability_bypass,
                diff_externals = true,
                npartials = n,
                static_size = true
            )
            push!(storage, s)
            b = InPlaceFlashBuffer(nc)
            push!(buffers, b)
        end
        storage = tuple(storage...)
        buffers = tuple(buffers...)
        new{typeof(method), typeof(storage), typeof(buffers)}(
            storage,
            method,
            buffers,
            threads,
            tolerance,
            tolerance_bypass,
            stability_bypass,
            reuse_guess
        )
    end
end

default_value(model, ::FlashResults) = FlashedMixture2Phase(model.system.equation_of_state)

function initialize_variable_value(model, pvar::FlashResults, val::AbstractDict; need_value = false, T = Jutul.float_type(model.context))
    @assert need_value == false
    n = number_of_entities(model, pvar)
    v = default_value(model, pvar)
    T = typeof(v)
    V = Vector{T}()
    sizehint!(V, n)
    for i in 1:n
        push!(V, default_value(model, pvar))
    end
    initialize_variable_value(model, pvar, V)
end

function Jutul.initialize_variable_value!(state, model, pvar::FlashResults, symb, val::AbstractDict; kwarg...)
    state[symb] = initialize_variable_value(model, pvar, val; kwarg...)
end

function initialize_variable_ad!(state, model, pvar::FlashResults, symb, npartials, diag_pos; context = DefaultContext(), kwarg...)
    n = number_of_entities(model, pvar)
    v_ad = get_ad_entity_scalar(1.0, npartials, diag_pos; kwarg...)
    ∂T = typeof(v_ad)
    eos = flow_system(model.system).equation_of_state

    r = FlashedMixture2Phase(eos, ∂T)
    T = typeof(r)
    V = Vector{T}(undef, n)
    for i in 1:n
        V[i] = FlashedMixture2Phase(eos, ∂T)
    end
    state[symb] = V
    return state
end

@jutul_secondary function update_flash!(flash_results, fr::FlashResults, model, Pressure, Temperature, OverallMoleFractions, ix)
    eos = model.system.equation_of_state
    flash_entity_loop!(flash_results, fr, model, eos, Pressure, Temperature, OverallMoleFractions, nothing, ix)
end

@jutul_secondary function update_flash!(flash_results, fr::FlashResults, model::LVCompositionalModel3Phase, Pressure, Temperature, OverallMoleFractions, ImmiscibleSaturation, ix)
    eos = model.system.equation_of_state
    flash_entity_loop!(flash_results, fr, model, eos, Pressure, Temperature, OverallMoleFractions, ImmiscibleSaturation, ix)
end

function flash_entity_loop!(flash_results, fr, model, eos, Pressure, Temperature, OverallMoleFractions, sw, ix)
    storage, m, buffers = fr.storage, fr.method, fr.update_buffer

    S, buf = thread_buffers(storage, buffers)
    update_flash_buffer!(buf, eos, Pressure, Temperature, OverallMoleFractions)
    buf_z = buf.z
    buf_forces = buf.forces
    flash_entity_loop_impl!(flash_results, ix, S, fr, m, eos, buf_z, buf_forces, Pressure, Temperature, OverallMoleFractions, sw)
end

function flash_entity_loop_impl!(flash_results, ix, S, fr, m, eos, buf_z, buf_forces, Pressure, Temperature, OverallMoleFractions, sw)
    extra_print = false
    if extra_print
        code_log = Dict(
            FLASH_SINGLE_PHASE_TESTED => 0,
            FLASH_SINGLE_PHASE_BYPASSED => 0,
            FLASH_FULL => 0,
            FLASH_RESTARTED => 0,
            FLASH_PURE_WATER => 0
        )
    end
    @inbounds for i in ix
        flash_results[i], code = internal_flash!(flash_results[i], S, fr, m, eos, buf_z, buf_forces, Pressure, Temperature, OverallMoleFractions, sw, i)
        if extra_print
            code_log[code] += 1
        end
    end
    if extra_print
        @info "Flash is done:" code_log
    end
end

@inline function Jutul.update_values!(vals::AbstractVector{<:FlashedMixture2Phase}, next::AbstractVector{<:FlashedMixture2Phase})
    function replace_flashed_phase_values(x, y)
        Jutul.update_values!(x.mole_fractions, y.mole_fractions)
        Z = x.Z
        return FlashedPhase(x.mole_fractions, Z - value(Z) + value(y.Z))
    end

    for i in eachindex(vals)
        old_v = vals[i]
        next_v = next[i]
        l = replace_flashed_phase_values(old_v.liquid, next_v.liquid)
        v = replace_flashed_phase_values(old_v.vapor, next_v.vapor)
        K = old_v.K
        V = old_v.V
        Jutul.update_values!(K, next_v.K)
        state = next_v.state
        V = V - value(V) + value(next_v.V)
        vals[i] = FlashedMixture2Phase(state, K, V, l, v)
    end
end

function thread_buffers(storage, buffers)
    thread_id = Threads.threadid()
    S = storage[thread_id]
    thread_buffer = buffers[thread_id]
    return (S, thread_buffer)
end

function update_flash_buffer!(buf, eos, Pressure, Temperature, OverallMoleFractions)
    maybe_ad_T = Base.promote_type(typeof(Pressure), typeof(Temperature), eltype(OverallMoleFractions))
    if isnothing(buf.forces) || eltype(buf.forces.A_ij) != maybe_ad_T
        P = Pressure[1]
        T = Temperature[1]
        if isa(OverallMoleFractions, AbstractVector)
            Z = OverallMoleFractions
        else
            Z = OverallMoleFractions[:, 1]
        end
        buf.forces = force_coefficients(eos, (p = P, T = T, z = Z), static_size = true)
    end
end

function update_flash_buffer!(buf, eos::KValuesEOS, Pressure, Temperature, OverallMoleFractions)
    return nothing
end

function internal_flash!(f, S, var, m, eos, buf_z, buf_forces, Pressure, Temperature, OverallMoleFractions, sw, i)
    @inline function immiscible_sat(::Nothing, i)
        return 0.0
    end

    @inline function immiscible_sat(s, i)
        return @inbounds s[i]
    end

    @inbounds begin
        P = Pressure[i]
        T = Temperature[i]
        Z = @view OverallMoleFractions[:, i]
        Sw = immiscible_sat(sw, i)

        K = f.K
        x = f.liquid.mole_fractions
        y = f.vapor.mole_fractions
        V = f.V
        b = f.critical_distance

        return update_flash_result(S, m, eos, f.state, K, f.flash_cond, f.flash_stability, x, y, buf_z, V, buf_forces, P, T, Z, Sw,
            critical_distance = b,
            tolerance = var.tolerance,
            stability_bypass = var.stability_bypass,
            reuse_guess = var.reuse_guess,
            tolerance_bypass = var.tolerance_bypass
        )
    end
end

import MultiComponentFlash: michelsen_critical_point_measure!
function update_flash_result(S, m, eos, phase_state, K, cond_prev, stability, x, y, z, V, forces, P, T, Z, Sw = 0.0;
        critical_distance::Float64 = NaN,
        tolerance::Float64 = 1e-8,
        tolerance_bypass::Float64 = 10.0,
        stability_bypass::Bool = false,
        reuse_guess::Bool = false,
    )

    # We can check for bypass if feature is enabled, we have a critical distance
    # that is finite and we were in single phase previously.
    # do_full_flash(c) = flash_2ph!(S, K, eos, c, NaN, method = m, extra_out = false, z_min = nothing)

    p_val = value(P)
    T_val = value(T)
    @. z = max(value(Z), MultiComponentFlash.MINIMUM_COMPOSITION)

    new_cond = (p = p_val, T = T_val, z = z)
    do_flash, critical_distance, single_phase_code = single_phase_bypass_check(eos, new_cond, cond_prev, phase_state, Sw, critical_distance, stability_bypass, tolerance_bypass)

    if do_flash
        vapor_frac, stability, return_code = two_phase_flash_implementation!(K, S, m, eos, phase_state, x, y, new_cond, V, reuse_guess)
        is_single_phase = isnan(vapor_frac)
        if is_single_phase && stability_bypass
            # If we did a full flash and single-phase prevailed outside the
            # shadow region it is time to update the critical distance for
            # future reference.
            if stability.liquid.trivial && stability.vapor.trivial
                # Outside shadow region, flash bypass can be enabled
                critical_distance = michelsen_critical_point_measure!(S.bypass, eos, new_cond.p, new_cond.T, new_cond.z)
            else
                # Inside shadow region, we cannot use flash bypass directly
                critical_distance = NaN
            end
        end
        # Update the condition only if we actually did a flash, otherwise we
        # keep the value where we last flashed around for stability bypass
        # testing.
        @. cond_prev.z = z
        cond_prev = (p = p_val, T = T_val, cond_prev.z)
    else
        return_code = single_phase_code
        is_single_phase = true
    end
    force_coefficients!(forces, eos, (p = P, T = T, z = Z))
    update_condition = false
    if is_single_phase
        # Single phase condition. Life is easy.
        Z_L, Z_V, V, phase_state = single_phase_update!(P, T, Z, x, y, forces, eos, new_cond)
    else
        # Two-phase condition: We have some work to do.
        Z_L, Z_V, V, phase_state = two_phase_update!(S, P, T, Z, x, y, K, vapor_frac, forces, eos, new_cond)
        # Reset critical distance to NaN and set condition to last flash since we are now in two-phase region.
        critical_distance = NaN
    end
    out = FlashedMixture2Phase(phase_state, K, V, x, y, Z_L, Z_V, critical_distance, cond_prev, stability)
    return (out, return_code)
end

function single_phase_bypass_check(eos, new_cond, old_cond, phase_state, Sw, critical_distance, stability_bypass, ϵ)
    is_pure_water = is_pure_single_phase(Sw)
    was_single_phase = phase_state == MultiComponentFlash.single_phase_l || phase_state == MultiComponentFlash.single_phase_v

    if is_pure_water
        # Set this to make sure that a proper stability test happens if we leave
        # the pure water condition.
        critical_distance = NaN
        need_two_phase_flash = false
        code = FLASH_PURE_WATER
    elseif stability_bypass && was_single_phase && isfinite(critical_distance)
        z_diff = p_diff = T_diff = 0.0
        # This is put early while we have the old z values in buffer
        for i in eachindex(old_cond.z, new_cond.z)
            z_diff = max(z_diff, abs(old_cond.z[i] - value(new_cond.z[i])))
        end
        p_diff = abs(new_cond.p - old_cond.p)
        T_diff = abs(new_cond.T - old_cond.T)

        b_old = critical_distance
        z_crit = z_diff ≥ b_old/ϵ
        p_crit = p_diff ≥ b_old*new_cond.p/ϵ
        T_crit = T_diff ≥ b_old*ϵ
        need_two_phase_flash = z_crit || p_crit || T_crit
        if need_two_phase_flash
            code = FLASH_SINGLE_PHASE_TESTED
        else
            code = FLASH_SINGLE_PHASE_BYPASSED
        end
        if false && need_two_phase_flash
            # Hacked in expensive test to check if the stability bypass is ok.
            # Can be manually activated.
            V = flash_2ph(eos, new_cond)
            if isnan(V)
                println("Bypass OK.")
            else
                error("Flash bypass was wrong for conditions $new_cond giving $V")
            end
        end
    else
        need_two_phase_flash = true
        code = FLASH_SINGLE_PHASE_TESTED
    end
    return (need_two_phase_flash, critical_distance, code)
end

function two_phase_flash_implementation!(K, S, m, eos, old_phase_state, x, y, flash_cond, V, reuse_guess)
    # Have to do some kind of flash, could be single or two-phase.
    need_full_flash = true
    if reuse_guess && old_phase_state == MultiComponentFlash.two_phase_lv
        # If the model was previously two-phase and this option is enabled, we
        # can try reusing the previous solution as an initial guess. This
        # follows Rasmussen et al where a theta estimating phase partition is
        # used, adapted to the K-value initial guess used by our flash.
        V = value(V)
        z = flash_cond.z
        K = estimate_K_values_from_previous_flash!(K, V, x, y, z)
        try
            vapor_frac, K, stats = MultiComponentFlash.flash_2ph_impl!(S, K, eos, flash_cond, V,
                method = SSINewtonFlash(swap_iter = 2),
                maxiter = 20,
                z_min = nothing,
                check = false
            )
            ϵ = 1e-6
            # Check for trivial solutions (i.e. phases are equal) and close to
            # single phase solutions. These are not valid flash results.
            trivial = true
            for k in K
                trivial = trivial && (abs(k-1.0) < ϵ)
            end
            tended_to_single_phase = vapor_frac <= ϵ || vapor_frac >= 1.0 - ϵ
            need_full_flash = !stats.converged || tended_to_single_phase || !isfinite(V) || trivial
        catch e
            jutul_message("Flash", "Exception ocurred in flash: $(typeof(e)), falling back to SSI with stability test.", color = :red)
        end
    end

    if need_full_flash
        vapor_frac, K, stats = MultiComponentFlash.flash_2ph_impl!(S, K, eos, flash_cond, NaN,
            method = m,
            z_min = nothing
        )
        code = FLASH_FULL
    else
        code = FLASH_RESTARTED
    end
    return (vapor_frac, stats.stability, code)
end

function estimate_K_values_from_previous_flash!(K, V, x, y, z)
    ϵ = MultiComponentFlash.MINIMUM_COMPOSITION
    x_t = 0.0
    y_t = 0.0
    for i in eachindex(z, x, y, K)
        z_i = z[i]
        x_i = value(x[i])
        y_i = value(y[i])

        # x_i = max(x_i, ϵ)
        # y_i = max(y_i, ϵ)
        θ_i = V*y_i/(V*y_i + (1.0 - V)*x_i)
        # Use K as working buffer
        K[i] = θ_i
        # Sum up molar amounts
        x_t += (1.0 - θ_i)*z[i]
        y_t += θ_i*z[i]
    end

    for i in eachindex(K)
        # Normalize estimates to mole fractions and get K-values
        θ_i = K[i]
        liquid = (1.0 - θ_i)/x_t
        vapor = θ_i/y_t
        K[i] = vapor/liquid # y_i / x_i
    end
    return K
end

function get_compressibility_factor(forces, eos, P, T, Z, phase = :unknown)
    ∂cond = (p = P, T = T, z = Z)
    force_coefficients!(forces, eos, ∂cond)
    scalars = force_scalars(eos, ∂cond, forces)
    return mixture_compressibility_factor(eos, ∂cond, forces, scalars, phase)
end

@inline function single_phase_update!(P, T, Z, x, y, forces, eos, c)
    AD = Base.promote_type(eltype(Z), typeof(P), typeof(T))
    Z_L = get_compressibility_factor(forces, eos, P, T, Z)
    Z_V = Z_L
    @. x = Z
    @. y = Z
    V = single_phase_label(eos, c)
    if V > 0.5
        phase_state = MultiComponentFlash.single_phase_v
    else
        phase_state = MultiComponentFlash.single_phase_l
    end
    V = convert(AD, V)
    out = (Z_L::AD, Z_V::AD, V::AD, phase_state::PhaseState2Phase)
    return out
end

function two_phase_pre!(S, P, T, Z, x::AbstractVector{AD}, y::AbstractVector{AD}, vapor_frac, eos, c) where {AD <: ForwardDiff.Dual}
    V = convert(AD, vapor_frac)
    if P isa AD || T isa AD || eltype(Z)<:AD
        inverse_flash_update!(S, eos, c, vapor_frac)
        ∂c = (p = P, T = T, z = Z)
        V = set_partials_vapor_fraction(V, S, eos, ∂c)
        set_partials_phase_mole_fractions!(x, S, eos, ∂c, :liquid)
        set_partials_phase_mole_fractions!(y, S, eos, ∂c, :vapor)
    end
    return V
end

two_phase_pre!(S, P, T, Z, x, y, V, eos, c) = V


@inline function two_phase_update!(S, P, T, Z, x, y, K, vapor_frac, forces, eos, c)
    AD = Base.promote_type(typeof(P), eltype(Z), typeof(T), eltype(x), eltype(y))
    @. x = liquid_mole_fraction(Z, K, vapor_frac)
    @. y = vapor_mole_fraction(x, K)
    V = two_phase_pre!(S, P, T, Z, x, y, vapor_frac, eos, c)
    Z_L = get_compressibility_factor(forces, eos, P, T, x, :liquid)
    Z_V = get_compressibility_factor(forces, eos, P, T, y, :vapor)
    phase_state = MultiComponentFlash.two_phase_lv

    return (convert(AD, Z_L), convert(AD, Z_V), convert(AD, V), phase_state)
end

function flash_entity_loop!(flash_results, fr, model, eos::KValuesEOS, Pressure, Temperature, OverallMoleFractions, sw, ix)
    storage, m, buffers = fr.storage, fr.method, fr.update_buffer
    S, buf = thread_buffers(storage, buffers)
    z_buf = buf.z
    kvalue_loop!(flash_results, ix, Pressure, Temperature, OverallMoleFractions, eos, z_buf)
end

function kvalue_loop!(flash_results, ix, Pressure, Temperature, OverallMoleFractions, eos, z_buf)
    @inbounds for i in ix
        P = Pressure[i]
        T = Temperature[i]
        Z = @view OverallMoleFractions[:, i]

        fr = flash_results[i]
        flash_results[i] = k_value_flash!(fr, eos, P, T, Z, z_buf)
    end
end

function k_value_flash!(result::FR, eos, P, T, Z, z) where FR
    Num_t = Base.promote_type(typeof(P), typeof(T), eltype(Z))
    ncomp = MultiComponentFlash.number_of_components(eos)
    c_ad = (p = P, T = T, z = Z)
    K_ad = eos.K_values_evaluator(c_ad)
    x = result.liquid.mole_fractions
    y = result.vapor.mole_fractions
    analytical_rr = ncomp == 2 || ncomp == 3
    K = result.K
    if analytical_rr
        # If we have 2 or 3 components the Rachford-Rice equations have an
        # analytical solution. We can then bypass a bunch of chain rule magic.
        V = flash_2ph!(nothing, K_ad, eos, c_ad, analytical = true)
    else
        @. z = max(value(Z), 1e-8)
        # Conditions
        c_numeric = (p = value(P), T = value(T), z = z)
        @inbounds for i in 1:ncomp
            K[i] = value(K_ad[i])
        end
        V = flash_2ph!(nothing, K, eos, c_numeric)
    end

    pure_liquid = V <= 0.0
    pure_vapor = V >= 1.0
    if pure_vapor || pure_liquid
        if pure_vapor
            phase_state = MultiComponentFlash.single_phase_v
        else
            phase_state = MultiComponentFlash.single_phase_l
        end
        @inbounds for i in 1:ncomp
            Z_i = Z[i]
            x[i] = Z_i
            y[i] = Z_i
        end
        V = Num_t(pure_vapor)
    else
        phase_state = MultiComponentFlash.two_phase_lv
        if !analytical_rr
            V = add_derivatives_to_vapor_fraction_rachford_rice(value(V), K_ad, Z, K, z)
        end
        V::Num_t
        @inbounds for i in 1:ncomp
            K_i = K_ad[i]
            x_i = liquid_mole_fraction(Z[i], K_i, V)
            x[i] = x_i
            y[i] = vapor_mole_fraction(x_i, K_i)
        end
    end
    V::Num_t
    Z_L = Z_V = convert(Num_t, 1.0)
    new_result = FlashedMixture2Phase(phase_state, K, V, x, y, Z_L, Z_V, NaN, result.flash_cond)
    return new_result::FR
end

function add_derivatives_to_vapor_fraction_rachford_rice(V::Float64, K, z, K_val = value(K), z_val = value(z))
    Zt = eltype(z)
    Kt = eltype(K)
    T = Base.promote_type(Zt, Kt)
    if T != Float64
        N = length(z)
        V0 = V
        V = convert(T, V)
        ∂R_chain = V.partials
        if Kt == T
            for i in 1:N
                dK_i = MultiComponentFlash.objectiveRR_dK(V0, K_val, z_val, i)
                ∂R_chain += K[i].partials*dK_i
            end
        end
        if Zt == T
            for i in 1:N
                dz_i = MultiComponentFlash.objectiveRR_dz(V0, K_val, z_val, i)
                ∂R_chain += z[i].partials*dz_i
            end
        end
        ∂R_dV = MultiComponentFlash.objectiveRR_dV(V0, K_val, z_val)
        ∂V = -∂R_chain/∂R_dV
        V = T(V0, ∂V)
    end
    return V
end
