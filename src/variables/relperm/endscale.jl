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
end

function EndPointScalingCoefficients(phases::Symbol)
    return EndPointScalingCoefficients{phases}()
end

Jutul.degrees_of_freedom_per_entity(model, ::EndPointScalingCoefficients) = 4

function Jutul.default_values(model, scalers::EndPointScalingCoefficients{P}) where P
    k = Symbol("scaler_$(P)_drainage")
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
        kr_i = JutulDarcy.table_by_region(kr, reg)
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
end

function add_scaling_parameters!(model::SimulationModel)
    add_scaling_parameters!(model.parameters, model[:RelativePermeabilities])
end

function add_scaling_parameters!(param, kr::AbstractRelativePermeabilities)
    if endpoint_scaling_is_active(kr)
        ph = kr.phases
        has_phase(x) = occursin("$x", "$ph")
        if has_phase(:w)
            param[:RelPermScalingW] = EndPointScalingCoefficients(:w)
        end
        if has_phase(:wo)
            param[:RelPermScalingOW] = EndPointScalingCoefficients(:ow)
        end
        if has_phase(:og)
            param[:RelPermScalingOG] = EndPointScalingCoefficients(:og)
        end
        if has_phase(:g)
            param[:RelPermScalingG] = EndPointScalingCoefficients(:g)
        end
    end
    return param
end

function get_kr_scalers(kr::PhaseRelativePermeability)
    return (kr.connate, kr.critical, kr.s_max, kr.k_max)
end

function get_kr_scalers(scaler::AbstractMatrix, c)
    @inbounds L = scaler[1, c]
    @inbounds CR = scaler[2, c]
    @inbounds U = scaler[3, c]
    @inbounds KM = scaler[4, c]
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


function get_endpoint_scalers(state, scaling::NoKrScale, ph)
    return nothing
end

function get_endpoint_scalers(state, scaling::Union{TwoPointKrScale, ThreePointKrScale}, ::Val{:wog})
    scalers = (
        w = state.RelPermScalingW,
        ow = state.RelPermScalingOW,
        og = state.RelPermScalingOG,
        g = state.RelPermScalingG
    )
    return scalers
end

function get_endpoint_scalers(state, scaling::Union{TwoPointKrScale, ThreePointKrScale}, ::Val{:wo})
    scalers = (
        w = state.RelPermScalingW,
        ow = state.RelPermScalingOW
    )
    return scalers
end

function get_endpoint_scalers(state, scaling::Union{TwoPointKrScale, ThreePointKrScale}, ::Val{:og})
    scalers = (
        og = state.RelPermScalingOG,
        g = state.RelPermScalingG
        )
    return scalers
end

function get_endpoint_scalers(state, scaling::Union{TwoPointKrScale, ThreePointKrScale}, ::Val{:wg})
    scalers = (
        w = state.RelPermScalingW,
        g = state.RelPermScalingG
    )
    return scalers
end
