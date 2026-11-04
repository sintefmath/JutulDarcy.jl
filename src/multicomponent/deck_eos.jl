"Deck cubic coefficients, retaining the underlying EOS alpha function."
struct DeckCubicParameters{E, A, B} <: MultiComponentFlash.AbstractGeneralizedCubic
    type::E
    omega_a::A
    omega_b::B
end

MultiComponentFlash.static_coefficients(t::DeckCubicParameters) = MultiComponentFlash.static_coefficients(t.type)

function underlying_deck_eos(eos::GenericCubicEOS{<:DeckCubicParameters})
    return GenericCubicEOS(eos.type.type, eos.mixture, eos.m_1, eos.m_2,
        eos.ω_a, eos.ω_b, eos.volume_shift)
end

function MultiComponentFlash.weight_ai(eos::GenericCubicEOS{<:DeckCubicParameters}, cond, i)
    base = underlying_deck_eos(eos)
    return MultiComponentFlash.weight_ai(base, cond, i)*(eos.type.omega_a[i]/eos.ω_a)
end

function MultiComponentFlash.weight_bi(eos::GenericCubicEOS{<:DeckCubicParameters}, cond, i)
    base = underlying_deck_eos(eos)
    return MultiComponentFlash.weight_bi(base, cond, i)*(eos.type.omega_b[i]/eos.ω_b)
end

MultiComponentFlash.forces_per_phase(eos::GenericCubicEOS{<:DeckCubicParameters}) = MultiComponentFlash.forces_per_phase(underlying_deck_eos(eos))

function MultiComponentFlash.binary_interaction(eos::GenericCubicEOS{<:DeckCubicParameters}, i::Int, j::Int, cond)
    return MultiComponentFlash.binary_interaction(underlying_deck_eos(eos), i, j, cond)
end

function MultiComponentFlash.static_cubic_type(t::DeckCubicParameters; float_type = missing)
    C = ismissing(float_type) ? identity : float_type
    base = MultiComponentFlash.static_cubic_type(t.type; float_type = float_type)
    return DeckCubicParameters(base, map(C, t.omega_a), map(C, t.omega_b))
end

function deck_cubic_coefficients(type, props, n)
    if !any(key -> haskey(props, key), ("OMEGAA", "OMEGAB"))
        return type
    end
    a, b, = MultiComponentFlash.static_coefficients(type)
    function coefficients(key, default)
        values = haskey(props, key) ? first(props[key]) : fill(default, n)
        length(values) == n || throw(ArgumentError("$key requires one value per component."))
        # Defaulted entries and the simulator's negative sentinel use the EOS default.
        values = map(x -> isnan(x) || x < 0 ? default : x, values)
        all(x -> isfinite(x) && x > 0, values) || throw(ArgumentError("$key coefficients must be finite and positive."))
        return Tuple(values)
    end
    return DeckCubicParameters(type, coefficients("OMEGAA", a), coefficients("OMEGAB", b))
end
