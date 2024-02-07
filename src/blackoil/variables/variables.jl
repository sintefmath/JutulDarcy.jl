Base.@kwdef struct GasMassFraction{R} <: ScalarVariable
    dz_max::R = 0.2
    sat_chop::Bool = true
end

maximum_value(::GasMassFraction) = 1.0
minimum_value(::GasMassFraction) = 1e-12
absolute_increment_limit(s::GasMassFraction) = s.dz_max

function update_primary_variable!(state, p::GasMassFraction, state_symbol, model, dx)
    s = state[state_symbol]
    update_gas_fraction!(s, state, p, model, dx)
end

struct BlackOilPhaseState <: ScalarVariable

end

Jutul.default_value(model, ::BlackOilPhaseState) = OilAndGas
Jutul.initialize_secondary_variable_ad!(state, model, var::BlackOilPhaseState, arg...; kwarg...) = state

struct Rs{F, R} <: ScalarVariable
    rs_max::F
    regions::R
end

function Rs(rs_max::F; regions::R = nothing) where {F, R}
    return Rs{F, R}(rs_max, regions)
end

function Jutul.subvariable(v::Rs, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(v.regions, c)
    return Rs(v.rs_max, regions = regions)
end

function Jutul.line_plot_data(model::SimulationModel, ::Rs)
    # TODO: Generalize to multiple regions
    (; X, F) = first(model.system.rs_max)
    return JutulLinePlotData(X[2:end]./1e5, F[2:end], title = "Saturated gas-in-oil ratio", xlabel = "Pressure [bar]", ylabel = "Rs")
end

struct Rv{F, R} <: ScalarVariable
    rv_max::F
    regions::R
end

function Rv(rv_max::F; regions::R = nothing) where {F, R}
    return Rv{F, R}(rv_max, regions)
end

function Jutul.subvariable(v::Rv, map::FiniteVolumeGlobalMap)
    c = map.cells
    regions = Jutul.partition_variable_slice(v.regions, c)
    return Rv(v.rv_max, regions = regions)
end

function Jutul.line_plot_data(model::SimulationModel, ::Rv)
    (; X, F) = first(model.system.rv_max)
    return JutulLinePlotData(X[2:end]./1e5, F[2:end], title = "Saturated oil-in-gas ratio", xlabel = "Pressure [bar]", ylabel = "Rv")
end

"""
    BlackOilUnknown(dr_max = Inf, ds_max = Inf)

Variable defining the variable switching black-oil variable. The `ds_max`
argument limits the maximum saturation change over a single Newton iteration
when both a oileic and gaseous phase is present. The `dr_max` limits the maximum
change in the undersaturated variable, taken relative to the maximum value of
the undersaturated variable.

During simulation, this variable can take on the following interpretations: Gas
saturation, Rs or Rv, depending on the phase conditions and miscibility model.
"""
Base.@kwdef struct BlackOilUnknown{R} <: ScalarVariable
    dr_max::R = Inf
    ds_max::R = 0.2
end

struct BlackOilX{T}
    val::T
    phases_present::PresentPhasesBlackOil
    sat_close::Bool
    function BlackOilX(val::T, phases = OilAndGas, sat_close = false) where T<:Real
        return new{T}(val, phases, sat_close)
    end
end

function Base.zero(::Type{BlackOilX{R}}) where R
    return BlackOilX(zero(R))
end

function Base.zero(::BlackOilX{R}) where R
    return BlackOilX(zero(R))
end

function Jutul.scalarized_primary_variable_type(model, var::BlackOilUnknown)
    return BlackOilX{Float64}
end

function Jutul.scalarize_primary_variable(model, source_vec, var::BlackOilUnknown, index)
    x = source_vec[index]
    return BlackOilX(value(x.val), x.phases_present, x.sat_close)
end

function Jutul.descalarize_primary_variable!(dest_array, model, V::BlackOilX, var::BlackOilUnknown, index)
    V_old = dest_array[index]
    val_converted = replace_value(V_old.val, V.val)
    dest_array[index] = BlackOilX(val_converted, V.phases_present, V.sat_close)
end

"""
    BlackOilX(sys::BlackOilVariableSwitchingSystem, p; sw = 0.0, so = 0.0, sg = 0.0, rs = 0.0, rv = 0.0, region = 1)

High level initializer for the black oil unknown degree of freedom. Will try to fill in the gaps unless system
is really underspecified.
"""
function BlackOilX(sys::BlackOilVariableSwitchingSystem, p; sw = 0.0, so = 0.0, sg = 0.0, rs = 0.0, rv = 0.0, region = 1)
    @assert p > 0 "Pressure must be positive"
    F_rs = sys.rs_max[region]
    F_rv = sys.rv_max[region]
    if has_disgas(sys)
        if sg > 0
            rs = F_rs(p)
            @assert sg ≈ 1 - sw
        else
            so = 1 - sw
        end
    end
    if has_vapoil(sys)
        if so > 0
            rv = F_rv(p)
            @assert so ≈ 1 - sw
        else
            sg = 1 - sw
        end
    end
    @assert sw + so + sg ≈ 1 "Saturations must sum up to one"

    return blackoil_unknown_init(F_rs, F_rv, sw, so, sg, rs, rv, p)
end

Jutul.default_value(model, ::BlackOilUnknown) = (NaN, OilAndGas, false) # NaN, Oil+Gas, away from bubble point
function Jutul.initialize_primary_variable_ad!(state, model, pvar::BlackOilUnknown, symb, npartials; offset, kwarg...)
    pre = state[symb]
    sys = model.system
    disgas = has_disgas(sys)
    vapoil = has_vapoil(sys)
    for (i, v) in enumerate(pre)
        if v.phases_present == GasOnly
            @assert vapoil "Cell $i initialized as GasOnly, but system does not have vapoil enabled."
        end
        if v.phases_present == OilOnly
            @assert disgas "Cell $i initialized as OilOnly, but system does not have vapoil enabled."
        end
    end
    vals = map(x -> x.val, pre)
    ad_vals = allocate_array_ad(vals, diag_pos = offset + 1, context = model.context, npartials = npartials; kwarg...)
    state[symb] = map((v, x) -> BlackOilX(v, x.phases_present, false), ad_vals, pre)
    return state
end

function Jutul.increment_norm(dX, state, model, X, pvar::BlackOilUnknown)
    T = eltype(dX)
    rs_max_v = rs_sum_v = zero(T)
    sg_max_v = sg_sum_v = zero(T)
    rv_max_v = rv_sum_v = zero(T)

    @inbounds for i in 1:length(dX)
        ph = X[i].phases_present
        dx_abs = abs(dX[i])
        if ph == OilAndGas
            sg_max_v = max(sg_max_v, dx_abs)
            sg_sum_v += dx_abs
        elseif ph == OilOnly
            rs_max_v = max(rs_max_v, dx_abs)
            rs_sum_v += dx_abs
        else
            rv_max_v = max(rv_max_v, dx_abs)
            rv_sum_v += dx_abs
        end
    end
    return (
            sg_sum = sg_sum_v, sg_max = sg_max_v,
            rs_sum = rs_sum_v, rs_max = rs_max_v,
            rv_sum = rv_sum_v, rv_max = rv_max_v,
            )
end

include("shrinkage.jl")
include("viscosity.jl")
include("density.jl")
include("total_masses.jl")

struct SurfaceVolumeMobilities <: PhaseVariables end

@jutul_secondary function update_surface_mob!(b_mob, var::SurfaceVolumeMobilities, model,
                                                        ShrinkageFactors,
                                                        PhaseMobilities,
                                                        ix)
    # For blackoil, the main upwind term
    for i in ix
        @inbounds for ph in axes(b_mob, 1)
            b_mob[ph, i] = ShrinkageFactors[ph, i]*PhaseMobilities[ph, i]
        end
    end
end

include("zg.jl")
include("varswitch.jl")
