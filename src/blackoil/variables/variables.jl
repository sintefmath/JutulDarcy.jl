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

struct Rs <: ScalarVariable end
function Jutul.line_plot_data(model::SimulationModel, ::Rs)
    (; X, F) = model.system.rs_max
    return JutulLinePlotData(X[2:end]./1e5, F[2:end], title = "Saturated gas-in-oil ratio", xlabel = "Pressure [bar]", ylabel = "Rs")
end

struct Rv <: ScalarVariable end

function Jutul.line_plot_data(model::SimulationModel, ::Rv)
    (; X, F) = model.system.rv_max
    return JutulLinePlotData(X[2:end]./1e5, F[2:end], title = "Saturated oil-in-gas ratio", xlabel = "Pressure [bar]", ylabel = "Rv")
end

Base.@kwdef struct BlackOilUnknown{R} <: ScalarVariable
    dr_max::R = Inf
    ds_max::R = 0.2
end

export BlackOilX
struct BlackOilX{T}
    val::T
    phases_present::PresentPhasesBlackOil
    sat_close::Bool
    function BlackOilX(val::T, phases = OilAndGas, sat_close = false) where T<:Real
        return new{T}(val, phases, sat_close)
    end
end

"""
    BlackOilX(sys::BlackOilVariableSwitchingSystem, p; sw = 0.0, so = 0.0, sg = 0.0, rs = 0.0, rv = 0.0)

High level initializer for the black oil unknown degree of freedom. Will try to fill in the gaps unless system
is really underspecified.
"""
function BlackOilX(sys::BlackOilVariableSwitchingSystem, p; sw = 0.0, so = 0.0, sg = 0.0, rs = 0.0, rv = 0.0)
    @assert p > 0 "Pressure must be positive"
    F_rs = sys.rs_max
    F_rv = sys.rv_max
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

include("shrinkage.jl")
include("viscosity.jl")
include("density.jl")
include("total_masses.jl")

struct SurfaceVolumeMobilities <: PhaseVariables end

@jutul_secondary function update_surface_mob!(b_mob, var::SurfaceVolumeMobilities, model,
                                                        ShrinkageFactors,
                                                        PhaseViscosities,
                                                        RelativePermeabilities,
                                                        ix)
    # For blackoil, the main upwind term
    for i in ix
        @inbounds for ph in axes(b_mob, 1)
            b_mob[ph, i] = ShrinkageFactors[ph, i]*RelativePermeabilities[ph, i]/PhaseViscosities[ph, i]
        end
    end
end

include("zg.jl")
include("varswitch.jl")
