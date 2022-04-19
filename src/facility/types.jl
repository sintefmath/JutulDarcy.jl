export TotalMassVelocityMassFractionsFlow

abstract type FacilitySystem <: JutulSystem end
struct PredictionMode <: FacilitySystem end
struct HistoryMode <: FacilitySystem end

abstract type SurfaceFacilityDomain <: JutulDomain end
abstract type WellControllerDomain <: SurfaceFacilityDomain end
struct WellGroup <: WellControllerDomain
    well_symbols::Vector{Symbol}
end

struct Wells <: JutulUnit end
struct TotalSurfaceMassRate <: ScalarVariable end
abstract type WellTarget end
abstract type SurfaceVolumeTarget <: WellTarget end

# Basics
export BottomHolePressureTarget, TotalRateTarget, SinglePhaseRateTarget
# Phase mixtures
export SurfaceLiquidRateTarget, SurfaceOilRateTarget, SurfaceWaterRateTarget, SurfaceGasRateTarget

struct BottomHolePressureTarget <: WellTarget
    value::AbstractFloat
end

struct SinglePhaseRateTarget <: SurfaceVolumeTarget
    value::AbstractFloat
    phase::AbstractPhase
end

lumped_phases(t::SinglePhaseRateTarget) = (t.phase, )

"""
Liquid rate (reservoir: oil + water but not gas)
"""
struct SurfaceLiquidRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end

lumped_phases(::SurfaceLiquidRateTarget) = (AqueousPhase(), LiquidPhase())

"""
Oil rate target
"""
struct SurfaceOilRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end

lumped_phases(::SurfaceOilRateTarget) = (LiquidPhase(), )

"""
Gas rate target
"""
struct SurfaceGasRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end

lumped_phases(::SurfaceGasRateTarget) = (VaporPhase(), )

"""
Water rate target
"""
struct SurfaceWaterRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end

lumped_phases(::SurfaceWaterRateTarget) = (AqueousPhase(), )

"""
All rates at surface conditions
"""
struct TotalRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end

struct DisabledTarget <: WellTarget end
abstract type WellForce <: JutulForce end
abstract type WellControlForce <: WellForce end

default_limits(ctrl) = as_limit(ctrl.target)
as_limit(target) = NamedTuple([Pair(translate_target_to_symbol(target, shortname = true), target.value)])
as_limit(T::DisabledTarget) = nothing

struct DisabledControl{T} <: WellControlForce
    target::T
    function DisabledControl()
        t = DisabledTarget()
        new{DisabledTarget}(t)
    end
end

function replace_target(f::DisabledControl, target)
    target::DisabledTarget()
    return f
end

struct InjectorControl{T} <: WellControlForce
    target::T
    injection_mixture
    mixture_density
    function InjectorControl(target::T, mix; density = 1.0) where T<:WellTarget
        if isa(mix, Real)
            mix = [mix]
        end
        mix = vec(mix)
        @assert sum(mix) ≈ 1
        new{T}(target, mix, density)
    end
end
replace_target(f::InjectorControl, target) = InjectorControl(target, f.injection_mixture, density = f.mixture_density)
default_limits(f::InjectorControl{T}) where T<:BottomHolePressureTarget = merge((rate_lower = MIN_ACTIVE_WELL_RATE, ), as_limit(f.target))


struct ProducerControl{T} <: WellControlForce
    target::T
    function ProducerControl(target::T) where T<:WellTarget
        new{T}(target)
    end
end

default_limits(f::ProducerControl{T}) where T<:SurfaceVolumeTarget = merge((bhp = 101325.0,), as_limit(f.target)) # 1 atm
default_limits(f::ProducerControl{T}) where T<:BottomHolePressureTarget = merge((rate_lower = -MIN_ACTIVE_WELL_RATE,), as_limit(f.target))

function replace_target(f::ProducerControl, target)
    return ProducerControl(target)
end

struct WellGroupConfiguration
    operating_controls # Currently operating control
    requested_controls # The requested control (which may be different if limits are hit)
    limits             # Operating limits for the wells
    function WellGroupConfiguration(well_symbols, control = nothing, limits = nothing)
        if isnothing(control)
            control = Dict{Symbol, WellControlForce}()
            for s in well_symbols
                control[s] = DisabledControl()
            end
        end
        requested = deepcopy(control)
        if isnothing(limits)
            limits = Dict{Symbol, Any}()
            for s in well_symbols
                limits[s] = nothing
            end
        end
        new(control, requested, limits)
    end
end

operating_control(cfg::WellGroupConfiguration, well::Symbol) = cfg.operating_controls[well]
current_limits(cfg::WellGroupConfiguration, well::Symbol) = cfg.limits[well]

struct ControlEquationWell <: JutulEquation
    # Equation:
    #        q_t - target = 0
    #        p|top cell - target = 0
    # We need to store derivatives with respect to q_t (same entity) and the top cell (other entity)
    equation::JutulAutoDiffCache
    function ControlEquationWell(model, number_of_equations; kwarg...)
        nw = count_entities(model.domain, Wells())
        alloc = (entity) -> CompactAutoDiffCache(number_of_equations, nw, model, entity = entity; kwarg...)
        target_well = alloc(Wells())
        new(target_well)
    end
end

struct TotalMassVelocityMassFractionsFlow <: FlowType end

struct PerforationMask{V} <: JutulForce where V<:AbstractVector
    values::V
    function PerforationMask(v::T) where T<:AbstractVecOrMat
        return new{T}(copy(vec(v)))
    end
end

import Base.copy
Base.copy(m::PerforationMask) = PerforationMask(copy(m.values))

translate_target_to_symbol(t; shortname = false) = Symbol(t)
translate_target_to_symbol(t::BottomHolePressureTarget; shortname = false) = shortname ? :bhp : Symbol("Bottom hole pressure")
translate_target_to_symbol(t::TotalRateTarget; shortname = false) = shortname ? :rate : Symbol("Surface total rate")
translate_target_to_symbol(t::SurfaceWaterRateTarget; shortname = false) = shortname ? :wrat : Symbol("Surface water rate")
translate_target_to_symbol(t::SurfaceLiquidRateTarget; shortname = false) = shortname ? :lrat : Symbol("Surface liquid rate (water + oil)")
translate_target_to_symbol(t::SurfaceOilRateTarget; shortname = false) = shortname ? :orat : Symbol("Surface oil rate")
translate_target_to_symbol(t::SurfaceGasRateTarget; shortname = false) = shortname ? :grat : Symbol("Surface gas rate")
