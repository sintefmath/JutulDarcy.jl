abstract type AbstractKrScale end
struct NoKrScale <: AbstractKrScale end
struct TwoPointKrScale <: AbstractKrScale end
struct ThreePointKrScale <: AbstractKrScale end

struct ScaledPhaseRelativePermeability{T, N, S<:AbstractKrScale} <: AbstractPhaseRelativePermeability{T, N}
    scaling::S
    unscaled_kr::PhaseRelativePermeability{T, N}
    "Connate saturation"
    connate::N
    "The saturation at which rel. perm. becomes positive"
    critical::N # Used
    "Largest value of rel. perm."
    k_max::N # Used 
    "Maximum saturation at which rel. perm. is k_max"
    s_max::N
    residual::N
    residual_base::N
end

function ScaledPhaseRelativePermeability(
        kr::PhaseRelativePermeability{T, M},
        scaling::AbstractKrScale;
        connate,
        critical,
        s_max,
        k_max,
        residual,
        residual_base
    ) where {T, M}
    # Promote types
    N = promote_type(
        M,
        typeof(connate),
        typeof(critical),
        typeof(s_max),
        typeof(k_max),
        typeof(residual),
        typeof(residual_base)
    )
    kr_conv = PhaseRelativePermeability{T, N}(
        kr.k,
        kr.label,
        N(kr.connate),
        N(kr.critical),
        N(kr.s_max),
        N(kr.k_max),
        N(kr.input_s_max)
    )
    return ScaledPhaseRelativePermeability{T, N, typeof(scaling)}(
        scaling,
        kr_conv,
        N(connate),
        N(critical),
        N(k_max),
        N(s_max),
        N(residual),
        N(residual_base)
    )
end

function (kr::ScaledPhaseRelativePermeability{T, N, S})(Sat) where {T, N, S<:NoKrScale}
    return kr.unscaled_kr(Sat)
end

function (kr::ScaledPhaseRelativePermeability{T, N, S})(s) where {T, N, S<:TwoPointKrScale}
    kr0 = kr.unscaled_kr
    s_scale = two_point_saturation_scaling(s, kr0.critical, kr.critical, kr0.s_max, kr.s_max)
    return (kr.k_max/kr0.k_max)*kr0(s_scale)
end

function (kr::ScaledPhaseRelativePermeability{T, N, S})(s) where {T, N, S<:ThreePointKrScale}
    kr0 = kr.unscaled_kr
    s_scale = three_point_saturation_scaling(s, kr0.critical, kr.critical, kr0.s_max, kr.s_max, kr.residual_base, kr.residual)
    return (kr.k_max/kr0.k_max)*kr0(s_scale)
end

struct EndPointScalingCoefficients{phases} <: VectorVariables
    drainage::Bool
end

"""
    EndPointScalingCoefficients(phases::Symbol; drainage = true)

Parameters that represent end-point scaling for a relative permeability phase
(drainage or imbibition) as parametrized by four points (connate, critical,
saturation at which maximum kr first occurs and maximum kr).
"""
function EndPointScalingCoefficients(phases::Symbol; drainage = true)
    return EndPointScalingCoefficients{phases}(drainage)
end

Jutul.degrees_of_freedom_per_entity(model, ::EndPointScalingCoefficients) = 4

Jutul.minimum_value(::EndPointScalingCoefficients) = 0.0
Jutul.maximum_value(::EndPointScalingCoefficients) = 1.0

function Jutul.default_values(model, scalers::EndPointScalingCoefficients{P}) where P
    if scalers.drainage
        k = Symbol("scaler_$(P)_drainage")
    else
        k = Symbol("scaler_$(P)_imbibition")
    end
    data_domain = model.data_domain
    nc = number_of_cells(data_domain)
    relperm = Jutul.get_variable(model, :RelativePermeabilities)
    kr = relperm[P]
    n = degrees_of_freedom_per_entity(model, scalers)
    if haskey(data_domain, k)
        kscale_model = data_domain[k]
        T = eltype(kscale_model)
    else
        T = Float64
    end
    kscale = zeros(T, n, nc)
    for i in 1:nc
        reg = JutulDarcy.region(relperm.regions, i)
        if scalers.drainage
            kr_i = table_by_region(kr, reg)
        else
            kr_i = imbibition_table_by_region(kr, reg)
        end
        (; connate, critical, s_max, k_max) = kr_i
        kscale[1, i] = connate
        kscale[2, i] = critical
        kscale[3, i] = s_max
        kscale[4, i] = k_max
    end
    my_isfinite(x) = isfinite(x)
    my_isfinite(x::Jutul.ST.ADval) = isfinite(x.val)
    if haskey(data_domain, k)
        kscale_model = data_domain[k]
        for i in eachindex(kscale, kscale_model)
            override = kscale_model[i]
            if my_isfinite(override)
                kscale[i] = override
            end
        end
    end
    return kscale
end

function endpoint_scaling_model(::AbstractRelativePermeabilities)
    return NoKrScale()
end

function endpoint_scaling_is_active(x::AbstractRelativePermeabilities)
    return endpoint_scaling_is_active(endpoint_scaling_model(x))
end

function endpoint_scaling_is_active(x::AbstractKrScale)
    return true
end

function endpoint_scaling_is_active(x::NoKrScale)
    return false
end

function add_scaling_parameters!(model::MultiModel)
    add_scaling_parameters!(reservoir_model(model))
    return model
end

function add_scaling_parameters!(model::SimulationModel)
    add_scaling_parameters!(model.parameters, model[:RelativePermeabilities])
    return model
end

function add_scaling_parameters!(param, kr::AbstractRelativePermeabilities)
    if endpoint_scaling_is_active(kr)
        hyst = hysteresis_is_active(kr)
        ph = kr.phases
        has_phase(x) = occursin("$x", "$ph")
        if has_phase(:w)
            param[:RelPermScalingW] = EndPointScalingCoefficients(:w, drainage = true)
            if hyst
                param[:RelPermScalingWi] = EndPointScalingCoefficients(:w, drainage = false)
            end
        end
        if has_phase(:wo)
            param[:RelPermScalingOW] = EndPointScalingCoefficients(:ow, drainage = true)
            if hyst
                param[:RelPermScalingOWi] = EndPointScalingCoefficients(:ow, drainage = false)
            end
        end
        if has_phase(:og)
            param[:RelPermScalingOG] = EndPointScalingCoefficients(:og, drainage = true)
            if hyst
                param[:RelPermScalingOGi] = EndPointScalingCoefficients(:og, drainage = false)
            end
        end
        if has_phase(:g)
            param[:RelPermScalingG] = EndPointScalingCoefficients(:g, drainage = true)
            if hyst
                param[:RelPermScalingGi] = EndPointScalingCoefficients(:g, drainage = false)
            end
        end
    end
    return param
end

function get_kr_scalers(kr::PhaseRelativePermeability)
    return (kr.connate, kr.critical, kr.s_max, kr.k_max)
end

function get_kr_scalers(scaler::AbstractMatrix, c)
    ϵ = 1e-4
    @inbounds L = scaler[1, c]
    @inbounds CR = scaler[2, c]
    @inbounds U = scaler[3, c]
    @inbounds KM = scaler[4, c]
    U = max(U, 2*ϵ)
    CR = min(CR, U - ϵ)
    L = min(L, CR - ϵ)
    KM = max(KM, 2*ϵ)
    return (L, CR, U, KM)
end

function relperm_scaling(scaling, F, s::T, cr, CR, u, U, km, KM, r, R) where T<:Real
    if scaling isa ThreePointKrScale
        S = three_point_saturation_scaling(s, cr, CR, u, U, r, R)
    else
        S = two_point_saturation_scaling(s, cr, CR, u, U)
    end
    return (KM/km)*F(S)
end

function three_point_saturation_scaling(s::T, cr, CR, u, U, r, R) where T<:Real
    # @assert r >= cr
    # @assert R >= CR
    # @assert u >= r
    # @assert U >= R
    if s < CR
        S = zero(T)
    elseif s >= CR && s < R
        S = (s - CR)*(r-cr)/(R-CR) + cr
    elseif s >= R && s < U
        S = (s - R)*(u-r)/(U-R) + r
    else
        S = one(T)
    end
    return S
end

function two_point_saturation_scaling(s::T, cr, CR, u, U) where T<:Real
    if s < CR
        S = zero(T)
    elseif s >= CR && s < U
        S = (s - CR)*(u-cr)/(U-CR) + cr
    else
        S = one(T)
    end
    return S
end


function get_endpoint_scalers(state, scaling::NoKrScale, ph; drainage::Bool = true)
    return nothing
end

function get_endpoint_scalers(state, scaling::Union{TwoPointKrScale, ThreePointKrScale}, ::Val{:wog}; drainage::Bool = true)
    if drainage
        scalers = (
            w = state.RelPermScalingW,
            ow = state.RelPermScalingOW,
            og = state.RelPermScalingOG,
            g = state.RelPermScalingG
        )
    else
        scalers = (
            w = state.RelPermScalingWi,
            ow = state.RelPermScalingOWi,
            og = state.RelPermScalingOGi,
            g = state.RelPermScalingGi
        )
    end
    return scalers
end

function get_endpoint_scalers(state, scaling::Union{TwoPointKrScale, ThreePointKrScale}, ::Val{:wo}; drainage::Bool = true)
    if drainage
        scalers = (
            w = state.RelPermScalingW,
            ow = state.RelPermScalingOW
        )
    else
        scalers = (
            w = state.RelPermScalingWi,
            ow = state.RelPermScalingOWi
        )
    end
    return scalers
end

function get_endpoint_scalers(state, scaling::Union{TwoPointKrScale, ThreePointKrScale}, ::Val{:og}; drainage::Bool = true)
    if drainage
        scalers = (
            og = state.RelPermScalingOG,
            g = state.RelPermScalingG
            )
    else
        scalers = (
            og = state.RelPermScalingOGi,
            g = state.RelPermScalingGi
        )
    end
    return scalers
end

function get_endpoint_scalers(state, scaling::Union{TwoPointKrScale, ThreePointKrScale}, ::Val{:wg}; drainage::Bool = true)
    if drainage
        scalers = (
            w = state.RelPermScalingW,
            g = state.RelPermScalingG
        )
    else
        scalers = (
            w = state.RelPermScalingWi,
            g = state.RelPermScalingGi
        )
    end
    return scalers
end

function get_three_phase_scaled_relperms(scaling, krw, krow, krog, krg, swcon, scaler_w, scaler_ow, scaler_og, scaler_g, c)
    L_w, CR_w, U_w, KM_w = get_kr_scalers(scaler_w, c)
    L_ow, CR_ow, U_ow, KM_ow = get_kr_scalers(scaler_ow, c)
    L_og, CR_og, U_og, KM_og = get_kr_scalers(scaler_og, c)
    L_g, CR_g, U_g, KM_g = get_kr_scalers(scaler_g, c)

    _, cr_w, u_w, km_w = get_kr_scalers(krw)
    l_w = swcon

    l_ow, cr_ow, u_ow, km_ow = get_kr_scalers(krow)
    l_og, cr_og, u_og, km_og = get_kr_scalers(krog)
    l_g, cr_g, u_g, km_g = get_kr_scalers(krg)

    # Residual water
    R_w = 1.0 - CR_ow - L_g
    r_w = 1.0 - cr_ow - l_g
    # Residual gas
    R_g = 1.0 - CR_og - L_w
    r_g = 1.0 - cr_og - l_w
    # Residual oil
    R_ow = 1.0 - CR_w - L_g
    r_ow = 1.0 - cr_w - l_g
    R_og = 1.0 - CR_g - L_w
    r_og = 1.0 - cr_w - l_w
    # Maximum saturation oil (in persence of water)
    U_ow = 1.0 - L_w - L_g
    u_ow = 1.0 - l_w - l_g
    # Maximum saturation oil (in persence of gas)
    U_og = 1.0 - L_g - L_w
    u_og = 1.0 - l_g - l_w

    Krw = ScaledPhaseRelativePermeability(krw, scaling, connate = L_w, critical = CR_w, k_max = KM_w, s_max = U_w, residual = R_w, residual_base = r_w)
    Krg = ScaledPhaseRelativePermeability(krg, scaling, connate = L_g, critical = CR_g, k_max = KM_g, s_max = U_g, residual = R_g, residual_base = r_g)
    Krow = ScaledPhaseRelativePermeability(krow, scaling, connate = L_ow, critical = CR_ow, k_max = KM_ow, s_max = U_ow, residual = R_ow, residual_base = r_ow)
    Krog = ScaledPhaseRelativePermeability(krog, scaling, connate = L_og, critical = CR_og, k_max = KM_og, s_max = U_og, residual = R_og, residual_base = r_og)

    return (Krw, Krow, Krog, Krg)
end


function get_two_phase_scaled_relperms(scaling, krw, krn, scaler_w, scaler_n, c)
    L_n, CR_n, U_n, KM_n = get_kr_scalers(scaler_n, c)
    l_n, cr_n, u_n, km_n = get_kr_scalers(krn)

    L_w, CR_w, U_w, KM_w = get_kr_scalers(scaler_w, c)
    l_w, cr_w, u_w, km_w = get_kr_scalers(krw)
    l_w = max(l_w, zero(l_w))

    R_w = 1.0 - CR_n
    r_w = 1.0 - cr_n

    R_n = 1.0 - CR_w
    r_n = 1.0 - cr_w
    U_n = 1.0 - L_w
    u_n = 1.0 - l_w

    Krw = ScaledPhaseRelativePermeability(krw, scaling, connate = L_w, critical = CR_w, k_max = KM_w, s_max = U_w, residual = R_w, residual_base = r_w)
    Krn = ScaledPhaseRelativePermeability(krn, scaling, connate = L_n, critical = CR_n, k_max = KM_n, s_max = U_n, residual = R_n, residual_base = r_n)

    return (Krw, Krn)
end

