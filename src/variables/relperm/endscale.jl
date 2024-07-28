@enum KrScale NoKrScale TwoPointKrScale ThreePointKrScale

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

scaling_type(::AbstractRelativePermeabilities) = NoKrScale

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
    if scaling == ThreePointKrScale
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
