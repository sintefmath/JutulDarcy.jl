export TotalMassVelocityMassFractionsFlow

abstract type FacilitySystem <: JutulSystem end
struct PredictionMode <: FacilitySystem end
struct HistoryMode <: FacilitySystem end

abstract type SurfaceFacilityDomain <: JutulDomain end
abstract type WellControllerDomain <: SurfaceFacilityDomain end
struct WellGroup <: WellControllerDomain
    well_symbols::Vector{Symbol}
end

struct Wells <: JutulEntity end
struct TotalSurfaceMassRate <: ScalarVariable end
abstract type WellTarget end
abstract type SurfaceVolumeTarget <: WellTarget end

"""
Perforations are connections from well cells to reservoir vcells
"""
struct Perforations <: JutulEntity end

struct WellIndices <: ScalarVariable end

Jutul.associated_entity(::WellIndices) = Perforations()
function Jutul.default_values(model, ::WellIndices)
    w = model.domain.grid
    return copy(w.perforations.WI)
end

Base.show(io::IO, t::SurfaceVolumeTarget) = print(io, "$(typeof(t)) with value $(t.value) [m^3/s] for $(join([typeof(p) for p in lumped_phases(t)], ", "))")

# Basics
export BottomHolePressureTarget, TotalRateTarget, SinglePhaseRateTarget, DisabledTarget
# Phase mixtures
export SurfaceLiquidRateTarget, SurfaceOilRateTarget, SurfaceWaterRateTarget, SurfaceGasRateTarget

struct BottomHolePressureTarget <: WellTarget
    value::AbstractFloat
end

"""
    SinglePhaseRateTarget(q, phase)

Single-phase well target with value `q` specified for `phase`.

# Examples
```julia-repl
julia> SinglePhaseRateTarget(0.001, LiquidPhase())
SinglePhaseRateTarget of 0.001 [m^3/s] for LiquidPhase()
```

"""
struct SinglePhaseRateTarget <: SurfaceVolumeTarget
    value::AbstractFloat
    phase::AbstractPhase
end

lumped_phases(t::SinglePhaseRateTarget) = (t.phase, )
"""
    SurfaceLiquidRateTarget(q)

Well target of specified liquid rate with value `q` (liquid/oil and water, but not gas)
at surface conditions.
"""
struct SurfaceLiquidRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end

lumped_phases(::SurfaceLiquidRateTarget) = (AqueousPhase(), LiquidPhase())

"""
    SurfaceOilRateTarget(q)

Well target of specified oil rate with value `q` at surface conditions.
"""
struct SurfaceOilRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end

lumped_phases(::SurfaceOilRateTarget) = (LiquidPhase(), )

"""
    SurfaceGasRateTarget(q)

Well target of specified gas rate with value `q` at surface conditions.
"""
struct SurfaceGasRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end

lumped_phases(::SurfaceGasRateTarget) = (VaporPhase(), )

"""
    SurfaceWaterRateTarget(q)

Well target of specified water rate with value `q` at surface conditions.
"""
struct SurfaceWaterRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end

lumped_phases(::SurfaceWaterRateTarget) = (AqueousPhase(), )

"""
    TotalRateTarget(q)

Well target of specified total rate of all phases with value `q` at surface conditions.
"""
struct TotalRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
end
Base.show(io::IO, t::TotalRateTarget) = print(io, "TotalRateTarget with value $(t.value) [m^3/s]")

"""
    DisabledTarget(q)

Disabled target used when a well is under `DisabledControl()` only.
"""
struct DisabledTarget <: WellTarget end
abstract type WellForce <: JutulForce end
abstract type WellControlForce <: WellForce end

"""
    default_limits(ctrl)

Create reasonable default limits for well control `ctrl`, for example to avoid BHP injectors turning into producers.
"""

default_limits(ctrl) = as_limit(ctrl.target)
as_limit(target) = NamedTuple([Pair(translate_target_to_symbol(target, shortname = true), target.value)])
as_limit(T::DisabledTarget) = nothing

"""
    DisabledControl()

Control that disables a well. If a well is disabled, it is disconnected from the surface network and no flow occurs
between the well and the top side. Mass transfer can still occur inside the well, and between the well and the reservoir.

See also [`ProducerControl`](@ref), [`InjectorControl`](@ref).
"""
struct DisabledControl{T} <: WellControlForce
    target::T
    function DisabledControl()
        t = DisabledTarget()
        new{DisabledTarget}(t)
    end
end

"""
    replace_target(ctrl, new_target)

Create new well control using `ctrl` as a template that operates under `new_target`.
"""
function replace_target

end

function replace_target(f::DisabledControl, target)
    target::DisabledTarget()
    return f
end

"""
    InjectorControl(target, mix, [density])

Well control that specifies injection into the reservoir. `target` specifies the type of target and `mix` defines the
injection mass fractions for all species in the model during injection. 

For example, for a three-component system made up of CO2, H2O and H2, setting [0.1, 0.6, 0.3] would mean
that the injection stream would contain 1 part CO2, 6 parts H2O and 3 parts H2 by mass. For an immiscible
system (e.g. `LiquidPhase(), VaporPhase()`) the species corresponds to phases and [0.3, 0.7] would mean a
3 to 7 mixture of liquid and vapor by mass.

The density of the injected fluid at surface conditions is given by `density` which is defaulted to 1.0
if not given.

See also [`ProducerControl`](@ref), [`DisabledControl`](@ref).
"""
struct InjectorControl{T} <: WellControlForce
    target::T
    injection_mixture
    mixture_density
    phases
    function InjectorControl(target::T, mix; density = 1.0, phases = ((1, 1.0),)) where T<:WellTarget
        if isa(mix, Real)
            mix = [mix]
        end
        mix = vec(mix)
        @assert sum(mix) ≈ 1
        new{T}(target, mix, density, phases)
    end
end
replace_target(f::InjectorControl, target) = InjectorControl(target, f.injection_mixture, density = f.mixture_density, phases = f.phases)
default_limits(f::InjectorControl{T}) where T<:BottomHolePressureTarget = merge((rate_lower = MIN_ACTIVE_WELL_RATE, ), as_limit(f.target))

"""
    ProducerControl(target)

Well control for production out of the reservoir. `target` specifies the type of target (for example `BottomHolePressureTarget()`).

See also [`DisabledControl`](@ref), [`InjectorControl`](@ref).
"""
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

function Jutul.update_values!(old::WellGroupConfiguration, new::WellGroupConfiguration)
    return WellGroupConfiguration(copy(new.operating_controls), copy(new.requested_controls), copy(new.limits))
end

operating_control(cfg::WellGroupConfiguration, well::Symbol) = cfg.operating_controls[well]
current_limits(cfg::WellGroupConfiguration, well::Symbol) = cfg.limits[well]

struct ControlEquationWell <: JutulEquation
    # Equation:
    #        q_t - target = 0
    #        p|top cell - target = 0
end

struct WellSegmentFlow{C, T<:AbstractVector} <: Jutul.FlowDiscretization
    cell_discretizations::C
    face_discretizations::T
    function WellSegmentFlow(well, z)
        # Face part
        N = get_neighborship(well)
        nf = size(N, 2)
        nc = number_of_cells(well)
        function F(i)
            l = N[1, i]
            r = N[2, i]
            gdz =  -gravity_constant*(z[l] - z[r])
            return (left = l, right = r, gdz = gdz, face = i)
        end
        fdisc = map(F, 1:nf)

        # Handle cell part
        cdisc = Jutul.half_face_map(N, nc)
        return new{typeof(cdisc), typeof(fdisc)}(cdisc, fdisc)
    end
end

Base.show(io::IO, t::MIME"text/plain", d::WellSegmentFlow) = print(io, "WellSegmentFlow")

function (D::WellSegmentFlow)(i, ::Faces)
    return D.face_discretizations[i]
end

function (D::WellSegmentFlow)(i, ::Cells)
    cd = D.cell_discretizations
    loc = cd.face_pos[i]:(cd.face_pos[i+1]-1)
    faces = @views cd.faces[loc]
    signs = @views cd.face_sign[loc]
    cells = @views cd.cells[loc]
    return (faces = cd.faces[loc], signs = signs, cells = cells)
end

export PerforationMask
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
