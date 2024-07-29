abstract type AbstractHysteresis end

Base.@kwdef struct KilloughHysteresis <: AbstractHysteresis
    tol::Float64 = 0.1
end

struct CarlsonHysteresis <: AbstractHysteresis end

struct JargonHysteresis <: AbstractHysteresis end

struct NoHysteresis <: AbstractHysteresis end

##
function kr_hysteresis(t, drain, imb, s, s_max)
    if s > s_max
        kr = drain(s)
    else
        kr = hysteresis_impl(t, drain, imb, s, s_max)
    end
    return kr
end

function hysteresis_impl(t::CarlsonHysteresis, drain, imb, s, s_max)
    kr_at_max = drain.k(s_max)
    s_meet = Jutul.linear_interp(imb.k.F, imb.k.X, kr_at_max)
    s_shifted = s + s_meet - s_max
    return imb.k(s_shifted)
end

function hysteresis_impl(h_model::KilloughHysteresis, drain, imb, S, S_max)
    S_crit_imbibition = imb.critical
    S_crit_drainage = drain.critical
    # TODO: Check that this matches that of imbibition?
    kr_s_max = drain.s_max
    K = 1.0/(S_crit_imbibition - S_crit_drainage) - 1.0/(kr_s_max - S_crit_drainage)
    M = 1.0 + h_model.tol*(kr_s_max - S_max)
    S_crit = S_crit_drainage + (S_max - S_crit_drainage)/(M + K*(S_max - S_crit_drainage))
    S_norm = S_crit_imbibition + (S - S_crit)*(kr_s_max - S_crit_imbibition)/(S_max - S_crit)
    kr = imb.k(S_norm)*drain.k(S_max)/drain.k(kr_s_max)
    return kr
end

function hysteresis_impl(h_model::JargonHysteresis, drain, imb, S, S_max)
    S_crit_imbibition = imb.critical
    S_crit_drainage = drain.critical
    # TODO: Check that this matches that of imbibition?
    kr_s_max = drain.s_max
    S_norm = S_crit_drainage + (S - S_crit_drainage)*(kr_s_max - S_crit_drainage)/(S_max - S_crit_drainage)
    kr = (imb.k(S_norm)/drain.k(S_norm))*drain.k(S)
    return kr
end
