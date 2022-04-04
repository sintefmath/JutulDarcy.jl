mutable struct InPlaceFlashBuffer
    z
    forces
    function InPlaceFlashBuffer(n)
        z = zeros(n)
        new(z, nothing)
    end
end

struct FlashResults <: ScalarVariable
    storage
    method
    update_buffer
    use_threads::Bool
    function FlashResults(system; method = SSIFlash(), threads = Threads.nthreads() > 1)
        eos = system.equation_of_state
        nc = MultiComponentFlash.number_of_components(eos)
        if has_other_phase(system)
            n = nc + 1
        else
            n = nc
        end
        storage = []
        buffers = []
        if threads
            N = Threads.nthreads()
        else
            N = 1
        end
        for i = 1:N
            s = flash_storage(eos, method = method, inc_jac = true, diff_externals = true, npartials = n, static_size = true)
            push!(storage, s)
            b = InPlaceFlashBuffer(nc)
            push!(buffers, b)
        end
        storage = tuple(storage...)
        buffers = tuple(buffers...)
        new(storage, method, buffers, threads)
    end
end

default_value(model, ::FlashResults) = FlashedMixture2Phase(model.system.equation_of_state)

function initialize_variable_value!(state, model, pvar::FlashResults, symb, val::AbstractDict; need_value = false)
    @assert need_value == false
    n = number_of_entities(model, pvar)
    v = default_value(model, pvar)
    T = typeof(v)
    V = Vector{T}()
    sizehint!(V, n)
    for i in 1:n
        push!(V, default_value(model, pvar))
    end
    initialize_variable_value!(state, model, pvar, symb, V)
end

function initialize_variable_ad(state, model, pvar::FlashResults, symb, npartials, diag_pos; context = DefaultContext(), kwarg...)
    n = number_of_entities(model, pvar)
    v_ad = get_ad_entity_scalar(1.0, npartials, diag_pos; kwarg...)
    ∂T = typeof(v_ad)
    eos = model.system.equation_of_state

    r = FlashedMixture2Phase(eos, ∂T)
    T = typeof(r)
    # T = MultiComponentFlash.flashed_mixture_array_type(eos, ∂T)
    V = Vector{T}(undef, n)
    for i in 1:n
        V[i] = FlashedMixture2Phase(eos, ∂T)
    end
    state[symb] = V
    return state
end

@jutul_secondary function update_as_secondary!(flash_results, fr::FlashResults, model, param, Pressure, Temperature, OverallMoleFractions)
    storage, m, buffers = fr.storage, fr.method, fr.update_buffer
    eos = model.system.equation_of_state

    for buf in buffers
        update_flash_buffer!(buf, eos, Pressure, Temperature, OverallMoleFractions)
    end
    perform_flash_for_all_cells!(flash_results, storage, m, eos, buffers, Pressure, Temperature, OverallMoleFractions, threads = fr.use_threads)
end


function perform_flash_for_all_cells!(flash_results, storage, m, eos, buffers, P, T, z; threads = true)
    flash_cell(i, S, buf) = internal_flash!(flash_results, S, m, eos, buf, P, T, z, i)
    if threads
        N = Threads.nthreads()
        @assert length(storage) == N == length(buffers)
        @inbounds @batch threadlocal=thread_buffers(storage, buffers) for i in eachindex(flash_results)
            # Unpack thread specific storage
            S, thread_buffer = threadlocal
            # Do flash
            flash_cell(i, S, thread_buffer)
        end
    else
        S, buf = thread_buffers(storage, buffers)
        @inbounds for i in eachindex(flash_results)
            flash_cell(i, S, buf)
        end
    end
end

function thread_buffers(storage, buffers)
    thread_id = Threads.threadid()
    S = storage[thread_id]
    thread_buffer = buffers[thread_id]
    return (S, thread_buffer)
end

function update_flash_buffer!(buf, eos, Pressure, Temperature, OverallMoleFractions)
    if isnothing(buf.forces) || eltype(buf.forces.A_ij) != eltype(OverallMoleFractions)
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

function internal_flash!(flash_results, S, m, eos, buf, Pressure, Temperature, OverallMoleFractions, i)
    # Ready to work
    @inbounds begin
        f = flash_results[i]
        P = Pressure[i]
        T = Temperature[i]
        Z = @view OverallMoleFractions[:, i]

        K = f.K
        x = f.liquid.mole_fractions
        y = f.vapor.mole_fractions

        flash_results[i] = update_flash_result(S, m, buf, eos, K, x, y, buf.z, buf.forces, P, T, Z)
    end
end


function update_flash_result(S, m, buffer, eos, K, x, y, z, forces, P, T, Z)
    @. z = max(value(Z), 1e-8)
    # Conditions
    c = (p = value(P), T = value(T), z = z)
    # Perform flash
    vapor_frac = flash_2ph!(S, K, eos, c, NaN, method = m, extra_out = false, z_min = nothing)
    force_coefficients!(forces, eos, (p = P, T = T, z = Z))
    if isnan(vapor_frac)
        # Single phase condition. Life is easy.
        Z_L, Z_V, V, phase_state = single_phase_update!(P, T, Z, x, y, forces, eos, c)
    else
        # Two-phase condition: We have some work to do.
        Z_L, Z_V, V, phase_state = two_phase_update!(S, P, T, Z, x, y, K, vapor_frac, forces, eos, c)
    end
    out = FlashedMixture2Phase(phase_state, K, V, x, y, Z_L, Z_V)
    return out
end

function get_compressibility_factor(forces, eos, P, T, Z)
    ∂cond = (p = P, T = T, z = Z)
    force_coefficients!(forces, eos, ∂cond)
    return mixture_compressibility_factor(eos, ∂cond, forces)
end

@inline function single_phase_update!(P, T, Z, x, y, forces, eos, c)
    AD = Base.promote_type(eltype(Z), typeof(P), typeof(T))
    Z_L = get_compressibility_factor(forces, eos, P, T, Z)
    Z_V = Z_L
    @. x = Z
    @. y = Z
    V = single_phase_label(eos.mixture, c)
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
    inverse_flash_update!(S, eos, c, vapor_frac)
    ∂c = (p = P, T = T, z = Z)
    V = set_partials_vapor_fraction(convert(AD, vapor_frac), S, eos, ∂c)
    set_partials_phase_mole_fractions!(x, S, eos, ∂c, :liquid)
    set_partials_phase_mole_fractions!(y, S, eos, ∂c, :vapor)
    return V
end

two_phase_pre!(S, P, T, Z, x, y, V, eos, c) = V


@inline function two_phase_update!(S, P, T, Z, x, y, K, vapor_frac, forces, eos, c)
    AD = Base.promote_type(typeof(P), eltype(Z), typeof(T))
    @. x = liquid_mole_fraction(Z, K, vapor_frac)
    @. y = vapor_mole_fraction(x, K)
    V = two_phase_pre!(S, P, T, Z, x, y, vapor_frac, eos, c)
    Z_L = get_compressibility_factor(forces, eos, P, T, x)
    Z_V = get_compressibility_factor(forces, eos, P, T, y)
    phase_state = MultiComponentFlash.two_phase_lv

    return (Z_L::AD, Z_V::AD, V::AD, phase_state)
end
