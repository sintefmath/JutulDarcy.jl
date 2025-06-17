
abstract type FacilitySystem <: JutulSystem end
struct PredictionMode <: FacilitySystem end
struct HistoryMode <: FacilitySystem end

const FacilityModel = SimulationModel{<:Any, <:FacilitySystem, <:Any, <:Any}

abstract type SurfaceFacilityDomain <: JutulDomain end
abstract type WellControllerDomain <: SurfaceFacilityDomain end
mutable struct WellGroup <: WellControllerDomain
    const well_symbols::Vector{Symbol} # Controlled wells
    "Can temporarily shut producers that try to reach zero rate multiple solves in a row"
    can_shut_producers::Bool
    "Can temporarily shut injectors that try to reach zero rate multiple solves in a row"
    can_shut_injectors::Bool
end

"""
    WellGroup(wells::Vector{Symbol}; can_shut_wells = true)

Create a well group that can control the given set of wells.
"""
function WellGroup(wells::Vector{Symbol}; can_shut_wells = true, can_shut_injectors = can_shut_wells, can_shut_producers = can_shut_wells)
    return WellGroup(wells, can_shut_producers, can_shut_injectors)
end

const WellGroupModel = SimulationModel{WellGroup, <:Any, <:Any, <:Any}

struct Wells <: JutulEntity end

"""
    TotalSurfaceMassRate(max_absolute_change = nothing, max_relative_change = nothing)

Variable, typically representing the primary variable for a [`WellGroup`](@ref).
The variable is a single entry per well and solves for the total surface mass
rate from a well to the facility model.
"""
Base.@kwdef struct TotalSurfaceMassRate <: ScalarVariable
    "Maximum absolute change betweeen two Newton updates (nominally kg/s)"
    max_absolute_change::Union{Float64, Nothing} = nothing
    "Maximum relative change between two Newton updates. Warning: Can be dangerous if set for wells operating around zero rate."
    max_relative_change::Union{Float64, Nothing} = nothing
end

Base.@kwdef struct SurfaceTemperature <: ScalarVariable
    "Maximum absolute change betweeen two Newton updates (nominally K)"
    max_absolute_change::Union{Float64, Nothing} = nothing
    "Maximum relative change between two Newton updates."
    max_relative_change::Union{Float64, Nothing} = nothing
    min = 273.15
    max = 1e6
end

function Jutul.absolute_increment_limit(q::TotalSurfaceMassRate)
    return q.max_absolute_change
end

function Jutul.relative_increment_limit(q::TotalSurfaceMassRate)
    return q.max_relative_change
end

abstract type WellTarget end
abstract type SurfaceVolumeTarget <: WellTarget end

"""
    Perforations()

Entity that defines perforations: Connections from well cells to reservoir
cells.
"""
struct Perforations <: JutulEntity end

"""
    WellIndices()

Parameter for the connection strength between a well and the reservoir for a
given perforation. Typical values come from a combination of Peaceman's formula,
upscaling and/or history matching.
"""
struct WellIndices <: ScalarVariable end

Jutul.minimum_value(::WellIndices) = 0.0
Jutul.variable_scale(::WellIndices) = 1e-10

Jutul.associated_entity(::WellIndices) = Perforations()
function Jutul.default_values(model, ::WellIndices)
    w = physical_representation(model.domain)
    return vec(copy(w.perforations.WI))
end

"""
    PerforationGravityDifference()

Parameter for the height difference from the wellbore and the connected node in
the well.
"""
struct PerforationGravityDifference <: ScalarVariable end

Jutul.associated_entity(::PerforationGravityDifference) = Perforations()
function Jutul.default_values(model, ::PerforationGravityDifference)
    w = physical_representation(model.domain)
    return vec(copy(w.perforations.gdz))
end

Base.show(io::IO, t::SurfaceVolumeTarget) = print(io, "$(typeof(t)) with value $(t.value) [m^3/s] for $(join([typeof(p) for p in lumped_phases(t)], ", "))")

"""
    BottomHolePressureTarget(q, phase)

Bottom-hole pressure (bhp) target with target pressure value `bhp`. A well
operating under a bhp constraint will keep the well pressure at the bottom hole
(typically the top of the perforations) fixed at this value unless doing so
would violate other constraints, like the well switching from injection to
production when declared as an injector.

# Examples
```julia-repl
julia> BottomHolePressureTarget(100e5)
BottomHolePressureTarget with value 100.0 [bar]
```

"""
struct BottomHolePressureTarget{T} <: WellTarget
    value::T
end
Base.show(io::IO, t::BottomHolePressureTarget) = print(io, "BottomHolePressureTarget with value $(convert_from_si(t.value, :bar)) [bar]")

"""
    SinglePhaseRateTarget(q, phase)

Single-phase well target with value `q` specified for `phase`.

# Examples
```julia-repl
julia> SinglePhaseRateTarget(0.001, LiquidPhase())
SinglePhaseRateTarget of 0.001 [m^3/s] for LiquidPhase()
```

"""
struct SinglePhaseRateTarget{T, P} <: SurfaceVolumeTarget
    value::T
    phase::P
end

lumped_phases(t::SinglePhaseRateTarget) = (t.phase, )
"""
    SurfaceLiquidRateTarget(q)

Well target of specified liquid rate at surface conditions with value `q`.
Typically used for a [`ProducerControl`](@ref) as you have full control over the
mixture composition during injection.

Liquid rate, sometimes abbreviated LRAT, is made up of the phases that remain
liquid at surface conditions. Typically, this will be water and oil if present
in the model, but never different types of gas. If a producing becomes nearly or
completely flooded by gas the well can go to very high or even infinite flows.
It is therefore important to combine this control with a limit such as a
bottom-hole-pressure constraint.
"""
struct SurfaceLiquidRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
    function SurfaceLiquidRateTarget(v::T) where T
        if T == Float64
            isfinite(v) || throw(ArgumentError("Rate must be finite, was $v"))
        end
        return new{T}(v)
    end
end

lumped_phases(::SurfaceLiquidRateTarget) = (AqueousPhase(), LiquidPhase())

"""
    SurfaceOilRateTarget(q)

Well target of specified oil rate with value `q` at surface conditions.
Typically used for a [`ProducerControl`](@ref) as oil, for economic reasons, is
rarely injected into the subsurface. Abbreviated as ORAT in some settings.
"""
struct SurfaceOilRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
    function SurfaceOilRateTarget(v::T) where T
        if T == Float64
            isfinite(v) || throw(ArgumentError("Rate must be finite, was $v"))
        end
        return new{T}(v)
    end
end

lumped_phases(::SurfaceOilRateTarget) = (LiquidPhase(), )

"""
    SurfaceGasRateTarget(q)

Well target of specified gas rate with value `q` at surface conditions.

Often used for both [`InjectorControl`](@ref) [`ProducerControl`](@ref).
Abbreviated as GRAT in some settings. If used for production it is important to
also impose limits, as the well rate may become very high if there is little gas
present.
"""
struct SurfaceGasRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
    function SurfaceGasRateTarget(v::T) where T
        if T == Float64
            isfinite(v) || throw(ArgumentError("Rate must be finite, was $v"))
        end
        return new{T}(v)
    end
end

lumped_phases(::SurfaceGasRateTarget) = (VaporPhase(), )

"""
    SurfaceWaterRateTarget(q)

Well target of specified water rate with value `q` at surface conditions.

Often used for both [`InjectorControl`](@ref) [`ProducerControl`](@ref). If used
for production it is important to also impose limits, as the well rate may
become very high if there is little water present.
"""
struct SurfaceWaterRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
    function SurfaceWaterRateTarget(v::T) where T
        if T == Float64
            isfinite(v) || throw(ArgumentError("Rate must be finite, was $v"))
        end
        return new{T}(v)
    end
end

lumped_phases(::SurfaceWaterRateTarget) = (AqueousPhase(), )

"""
    TotalRateTarget(q)

Well target of specified total rate (sum of all phases) with value `q` at surface
conditions.

Often used for both [`InjectorControl`](@ref) [`ProducerControl`](@ref).
"""
struct TotalRateTarget{T} <: SurfaceVolumeTarget where T<:AbstractFloat
    value::T
    function TotalRateTarget(v::T) where T
        if T == Float64
            isfinite(v) || throw(ArgumentError("Rate must be finite, was $v"))
        end
        return new{T}(v)
    end
end
Base.show(io::IO, t::TotalRateTarget) = print(io, "TotalRateTarget with value $(t.value) [m^3/s]")

"""
    TotalReservoirRateTarget(q)

Well target of specified total rate (sum of all phases) with value `q` at reservoir
conditions.

Often used for both [`InjectorControl`](@ref) [`ProducerControl`](@ref).
"""
struct TotalReservoirRateTarget{T} <: WellTarget where T<:AbstractFloat
    value::T
    function TotalReservoirRateTarget(v::T) where T
        if T == Float64
            isfinite(v) || throw(ArgumentError("Rate must be finite, was $v"))
        end
        return new{T}(v)
    end
end
Base.show(io::IO, t::TotalReservoirRateTarget) = print(io, "TotalReservoirRateTarget with value $(t.value) [m^3/s]")

"""
    HistoricalReservoirVoidageTarget(q, weights)

Historical RESV target for history matching cases. See
[`ReservoirVoidageTarget`](@ref). For historical rates, the weights described in
that target are computed based on the reservoir pressure and conditions at the
previous time-step.
"""
struct HistoricalReservoirVoidageTarget{T, K} <: WellTarget where {T<:AbstractFloat, K<:Tuple}
    value::T
    weights::K
end
Base.show(io::IO, t::HistoricalReservoirVoidageTarget) = print(io, "HistoricalReservoirVoidageTarget with value $(t.value) [m^3/s]")

"""
    ReservoirVoidageTarget(q, weights)

RESV target for history matching cases. The `weights` input should
have one entry per phase (or pseudocomponent) in the system. The well control
equation is then:

``|q_{ctrl} - \\sum_i w_i q_i^s|``

where ``q_i^s`` is the surface rate of phase ``i`` and ``w_i`` the weight of
component stream ``i``.

This constraint is typically set up from .DATA files for black-oil and immiscible cases.
"""
struct ReservoirVoidageTarget{T, K} <: WellTarget where {T<:AbstractFloat, K<:Tuple}
    value::T
    weights::K
end

mutable struct ReinjectionTarget <: WellTarget
    value
    wells::Vector{Symbol}
end

"""
    DisabledTarget(q)

Disabled target used when a well is under `DisabledControl()` only. The well
will be disconnected from the surface.
"""
struct DisabledTarget <: WellTarget end
abstract type WellForce <: JutulForce end
abstract type WellControlForce <: WellForce end

"""
    default_limits(ctrl)

Create reasonable default limits for well control `ctrl`, for example to avoid
BHP injectors turning into producers.
"""
function default_limits(ctrl)
    as_limit(ctrl.target)
end

as_limit(target) = NamedTuple([Pair(translate_target_to_symbol(target, shortname = true), target.value)])
as_limit(T::DisabledTarget) = nothing
as_limit(T::HistoricalReservoirVoidageTarget) = nothing
as_limit(target::ReservoirVoidageTarget) = NamedTuple([Pair(translate_target_to_symbol(target, shortname = true), (target.value, target.weights))])
as_limit(target::ReinjectionTarget) = NamedTuple([Pair(translate_target_to_symbol(target, shortname = true), (Inf))])

"""
    DisabledControl()

Control that disables a well. If a well is disabled, it is disconnected from the
surface network and no flow occurs between the well and the top side. Mass
transfer can still occur inside the well, and between the well and the reservoir
unless perforations are also closed by a [`PerforationMask`](@ref).

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
    target::DisabledTarget
    return f
end

function update_target!(ctrl, target, state_facility, state_well, facility)
    nothing
end

function update_target!(ctrl, target::ReinjectionTarget, state_facility, state_well, facility)

    q = 0.0
    for w in target.wells
        pos = get_well_position(facility.domain, w)
        qw = state_facility.TotalSurfaceMassRate[pos]
        @assert qw <= 0.0
        q -= qw
    end
    ρ = ctrl.mixture_density

    value = max(q./ρ, 1e-10)
    target.value = value

end

"""
    InjectorControl(target, mix, density = 1.0, phases = ((1, 1.0)), temperature = 293.15)

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
struct InjectorControl{T, R, P, M, E, TR} <: WellControlForce
    target::T
    injection_mixture::M
    mixture_density::R
    phases::P
    temperature::R
    enthalpy::E
    factor::R
    tracers::TR
    function InjectorControl(target::T, mix;
            density::Real = 1.0,
            phases = ((1, 1.0),),
            temperature::Real = 293.15,
            enthalpy = missing,
            tracers = missing,
            check = true,
            factor::Real = 1.0
        ) where {T<:WellTarget}
        density, temperature, factor = promote(density, temperature, factor)
        R = typeof(density)
        if isa(mix, Real)
            mix = [mix]
        end
        mix = vec(mix)
        if check && R == Float64
            @assert sum(mix) ≈ 1
            @assert isfinite(density) && density > 0.0 "Injector density must be finite and positive"
            @assert isfinite(temperature) && temperature > 0.0 "Injector temperature must be finite and positive"
        end
        if isa(tracers, Real)
            tracers = [tracers]
        end
        new{T, R, typeof(phases), typeof(mix), typeof(enthalpy), typeof(tracers)}(target, mix, density, phases, temperature, enthalpy, factor, tracers)
    end
end

function replace_target(f::InjectorControl, target, temperature = f.temperature)
    _, fact, den, T = Base.promote(target.value, f.factor, f.mixture_density, temperature)
    return InjectorControl(
        target,
        f.injection_mixture,
        temperature = T,
        enthalpy = f.enthalpy,
        density = den,
        phases = f.phases,
        factor = fact,
        tracers = f.tracers,
        check = false
    )
end

default_limits(f::InjectorControl{T}) where T<:BottomHolePressureTarget = merge((rate_lower = MIN_ACTIVE_WELL_RATE, ), as_limit(f.target))

function Base.isequal(f::InjectorControl, g::InjectorControl)
    t_eq = f.target == g.target 
    mix_eq = f.injection_mixture == g.injection_mixture
    den_eq = f.mixture_density == g.mixture_density
    phases_eq = f.phases == g.phases
    t_eq = f.temperature == g.temperature
    f_eq = f.factor == g.factor
    if ismissing(f.enthalpy)
        e_eq = ismissing(g.enthalpy)
    else
        e_eq = f.enthalpy == g.enthalpy
    end
    if ismissing(f.tracers)
        tr_eq = ismissing(g.tracers)
    else
        tr_eq = f.tracers == g.tracers
    end
    return t_eq && mix_eq && den_eq && phases_eq && t_eq && f_eq && e_eq && tr_eq
end

"""
    ProducerControl(target)

Well control for production out of the reservoir. `target` specifies the type of target (for example `BottomHolePressureTarget()`).

See also [`DisabledControl`](@ref), [`InjectorControl`](@ref).
"""
struct ProducerControl{T, R} <: WellControlForce
    target::T
    factor::R
    function ProducerControl(target::T; factor::R = 1.0) where {T<:WellTarget, R<:Real}
        new{T, R}(target, factor)
    end
end

default_limits(f::ProducerControl{T}) where T<:SurfaceVolumeTarget = merge((bhp = DEFAULT_MINIMUM_PRESSURE,), as_limit(f.target)) # 1 atm
default_limits(f::ProducerControl{T}) where T<:BottomHolePressureTarget = merge((rate_lower = -MIN_ACTIVE_WELL_RATE,), as_limit(f.target))

function replace_target(f::ProducerControl, target)
    return ProducerControl(target, factor = f.factor)
end

effective_surface_rate(qts, ::DisabledControl) = qts
effective_surface_rate(qts, c::Union{InjectorControl, ProducerControl}) = qts*c.factor

mutable struct WellGroupConfiguration{T, O, L}
    const operating_controls::T # Currently operating control
    const requested_controls::O # The requested control (which may be different if limits are hit)
    const limits::L             # Operating limits for the wells
    step_index::Int             # Internal book-keeping of what step we are at
    function WellGroupConfiguration(; operating, limits, requested = operating, step = 0)
        new{typeof(operating), typeof(requested), typeof(limits)}(operating, requested, limits, step)
    end
end

function WellGroupConfiguration(well_symbols, control = nothing, limits = nothing, step = 0)
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
    return WellGroupConfiguration(
        operating = control,
        requested = requested,
        limits = limits,
        step = step
    )
end

function Base.copy(c::WellGroupConfiguration)
    return WellGroupConfiguration(
        operating = copy(c.operating_controls),
        requested = copy(c.requested_controls),
        limits = copy(c.limits),
        step = c.step_index
    )
end

function Jutul.numerical_type(tc::WellGroupConfiguration)
    return Float64
end

function Jutul.update_values!(old::WellGroupConfiguration, new::WellGroupConfiguration)
    for (k, v) in new.operating_controls
        old.operating_controls[k] = v
    end
    for (k, v) in new.requested_controls
        old.requested_controls[k] = v
    end
    for (k, v) in new.limits
        old.limits[k] = v
    end
    old.step_index = new.step_index
    return old
end

operating_control(cfg::WellGroupConfiguration, well::Symbol) = cfg.operating_controls[well]
current_limits(cfg::WellGroupConfiguration, well::Symbol) = cfg.limits[well]

struct ControlEquationWell <: JutulEquation
    # Equation:
    #        q_t - target = 0
    #        p|top cell - target = 0
end

struct SurfaceTemperatureEquation <:JutulEquation
    # Equation:
    #        T_surf - T|top_cell = 0
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

"""
    mask = PerforationMask(mask::Vector)

Create a perforation mask. This can be passed to [`setup_forces`](@ref) for a
well under the `mask` argument. The mask should equal the number of perforations
in the well and is applied to the reference well indices in a multiplicative
fashion. For example, if a well named `:Injector` has two perforations, the
following mask would disable the first perforation and decrease the connection
strength for the second perforation by 50%:
```julia
mask = PerforationMask([0.0, 0.5])
iforces = setup_forces(W, mask = mask)
forces = setup_reservoir_forces(model, control = controls, Injector = iforces)
```
"""
struct PerforationMask{V} <: JutulForce where V<:AbstractVector
    values::V
    function PerforationMask(v::T) where T<:AbstractVecOrMat
        return new{T}(copy(vec(v)))
    end
end

import Base.copy
Base.copy(m::PerforationMask) = PerforationMask(copy(m.values))

function translate_target_to_symbol(t::T; shortname = true) where T
    info = well_target_information(t)
    if ismissing(info)
        ret = Symbol(T)
    elseif shortname
        ret = info.symbol
    else
        ret = Symbol(info.description)
    end
    return ret::Symbol
end

function well_target_information(;
        symbol::Symbol,
        description::String,
        unit_type::Symbol,
        unit_label::String,
        explanation::String = description,
        is_rate::Bool = true
    )
    return (
        symbol = symbol,
        description = description,
        explanation = explanation,
        unit_label = unit_label,
        unit_type = unit_type,
        is_rate = is_rate
    )
end

function well_target_information(x)
    return missing
end

function well_target_information(x::Symbol)
    if endswith("$x", "_mass_rate")
        cname = "$x"[1:end-10]
        out = well_target_information(
            symbol = x,
            description = "Component mass rate for $cname component",
            unit_type = :mass,
            unit_label = "kg/s",
            is_rate = true
        )
    else
        well_target_information(Val(x))
    end
end

function well_target_information(t::Union{BottomHolePressureTarget, Val{:bhp}})
    return well_target_information(
        symbol = :bhp,
        description = "Bottom hole pressure",
        explanation = "Pressure at well bottom hole. This is often given at or near the top perforation, but can be manually set to other depths.",
        unit_type = :pressure,
        unit_label = "Pa",
        is_rate = false
    )
end

function well_target_information(t::Union{TotalRateTarget, Val{:rate}})
    return well_target_information(
        symbol = :rate,
        description = "Surface total rate",
        explanation = "Total volumetric rate at surface conditions. This is the sum of all phases. For most models, it is the sum of the mass rates divided by the prescribed surface densities. For compositional models the density is computed using a flash.",
        unit_type = :liquid_volume_surface,
        unit_label = "m³/s"
    )
end


function well_target_information(t::Union{TotalReservoirRateTarget, Val{:resv_rate}})
    return well_target_information(
        symbol = :resv_rate,
        description = "Surface total rate",
        explanation = "Total volumetric rate at reservoir conditions. This is the sum of all phases. For most models, it is the sum of the mass rates divided by the prescribed surface densities.",
        unit_type = :liquid_volume_reservoir,
        unit_label = "m³/s"
    )
end

function well_target_information(t::Union{SurfaceWaterRateTarget, Val{:wrat}})
    return well_target_information(
        symbol = :wrat,
        description = "Surface water rate",
        explanation = "Water volumetric rate at surface conditions. This is the water mass stream divided by the surface density of water, which is typically around 1000 kg/m³",
        unit_type = :liquid_volume_surface,
        unit_label = "m³/s"
    )
end

function well_target_information(t::Union{SurfaceLiquidRateTarget, Val{:lrat}})
    return well_target_information(
        symbol = :lrat,
        description = "Surface water rate",
        explanation = "Liquid volumetric rate at surface conditions. This is the sum of the oil rate and the water rate.",
        unit_type = :liquid_volume_surface,
        unit_label = "m³/s"
    )
end

function well_target_information(t::Union{SurfaceOilRateTarget, Val{:orat}})
    return well_target_information(
        symbol = :orat,
        description = "Surface oil rate",
        explanation = "Oil rate at surface conditions. This is oil mass rate divided by the surface density of the oil phase.",
        unit_type = :liquid_volume_surface,
        unit_label = "m³/s"
    )
end

function well_target_information(t::Union{SurfaceGasRateTarget, Val{:grat}})
    return well_target_information(
        symbol = :grat,
        description = "Surface gas rate",
        explanation = "Gas rate at surface conditions. This is gas mass rate divided by the surface density of the gas phase.",
        unit_type = :gas_volume_surface,
        unit_label = "m³/s"
    )
end

function well_target_information(t::Union{ReservoirVoidageTarget, Val{:resv}})
    return well_target_information(
        symbol = :resv,
        description = "Reservoir voidage rate",
        explanation = "Reservoir voidage rate corresponds to a rate given at averaged pressure at reservoir conditions.",
        unit_type = :liquid_volume_reservoir,
        unit_label = "m³/s"
    )
end

function well_target_information(t::Union{HistoricalReservoirVoidageTarget, Val{:resv_history}})
    return well_target_information(
        symbol = :resv_history,
        description = "Historical reservoir voidage rate",
        explanation = "Historical reservoir voidage rate is a special rate used to match observed production rates.",
        unit_type = :liquid_volume_reservoir,
        unit_label = "m³/s"
    )
end

function well_target_information(t::Val{:mass_rate})
    return well_target_information(
        symbol = :mass_rate,
        description = "Total mass rate",
        explanation = "Total mass rate passing into or out from reservoir from well.",
        unit_type = :mass,
        unit_label = "kg/s"
    )
end

function well_target_information(t::Val{:control})
    return well_target_information(
        symbol = :control,
        description = "Control",
        explanation = "Type of control in use by well.",
        unit_type = :none,
        unit_label = "-",
        is_rate = false
    )
end

function well_target_information(t::Val{:temperature})
    return well_target_information(
        symbol = :temperature,
        description = "Well temperature",
        explanation = "Temperature at well bottom hole. This is often given at or near the top perforation, but can be manually set to other depths.",
        unit_type = :absolute_temperature,
        unit_label = "°K",
        is_rate = false
    )
end

function well_target_information(t::Val{:gor})
    return well_target_information(
        symbol = :gor,
        description = "Gas-oil-ratio",
        explanation = "Gas-oil ratio of production stream at surface conditions",
        unit_type = :id,
        unit_label = "",
        is_rate = false
    )
end

function well_target_information(t::Val{:wcut})
    return well_target_information(
        symbol = :wcut,
        description = "Water cut",
        explanation = "Volume fraction water in liquid production stream at surface conditions",
        unit_type = :id,
        unit_label = "",
        is_rate = false
    )
end

function realize_control_for_reservoir(state, ctrl, model, dt)
    return (ctrl, false)
end

function realize_control_for_reservoir(rstate, ctrl::ProducerControl{<:HistoricalReservoirVoidageTarget}, model, dt)
    sys = model.system
    w = ctrl.target.weights
    pv_t = 0.0
    p_avg = 0.0
    rs_avg = 0.0
    rv_avg = 0.0
    disgas = has_disgas(sys)
    vapoil = has_vapoil(sys)
    for c in eachindex(rstate.Pressure)
        p = value(rstate.Pressure[c])
        vol = value(rstate.FluidVolume[c])
        sw = value(rstate.ImmiscibleSaturation[c])
        vol_hc = vol*(1.0 - sw)
        p_avg += p*vol_hc
        pv_t += vol_hc
        if disgas
            rs_avg += value(rstate.Rs[c])*vol_hc
        end
        if vapoil
            rv_avg += value(rstate.Rv[c])*vol_hc
        end
    end
    p_avg /= pv_t
    rs_avg /= pv_t
    rv_avg /= pv_t

    a, l, v = phase_indices(sys)
    ww = w[a]
    wo = w[l]
    wg = w[v]
    if wo <= 1e-20
        rs = 0.0
    else
        rs = min(wg/wo, rs_avg)
    end
    if wg <= 1e-20
        rv = 0.0
    else
        rv = min(wo/wg, rv_avg)
    end

    svar = Jutul.get_secondary_variables(model)
    b_var = svar[:ShrinkageFactors]
    reg = b_var.regions
    bW = shrinkage(b_var.pvt[a], reg, p_avg, 1)
    if has_disgas(model.system)
        bO = shrinkage(b_var.pvt[l], reg, p_avg, rs, 1)
    else
        bO = shrinkage(b_var.pvt[l], reg, p_avg, 1)
    end
    if has_vapoil(model.system)
        bG = shrinkage(b_var.pvt[v], reg, p_avg, rv, 1)
    else
        bG = shrinkage(b_var.pvt[v], reg, p_avg, 1)
    end

    shrink = max(1.0 - rs*rv, 1e-20)
    shrink_avg = max(1.0 - rs_avg*rv_avg, 1e-20)
    old_rate = ctrl.target.value
    # Water
    new_water_rate = old_rate*ww/bW
    new_water_weight = 1/bW
    # Oil
    qo = old_rate*wo
    new_oil_rate = qo
    new_oil_weight = 1.0/(bO*shrink_avg)
    # Gas
    qg = old_rate*wg
    new_gas_rate = qg
    new_gas_weight = 1.0/(bG*shrink_avg)
    # Miscibility adjustments
    if vapoil
        new_oil_rate -= rv*qg
        new_gas_weight -= rv/(bO*shrink_avg)
    end
    if disgas
        new_gas_rate -= rs*qo
        new_oil_weight -= rs/(bG*shrink_avg)
    end
    new_oil_rate /= (bO*shrink)
    new_gas_rate /= (bG*shrink)


    new_weights = (new_water_weight, new_oil_weight, new_gas_weight)
    @assert all(isfinite, new_weights) "Computed RESV weights were non-finite: $new_weights"

    new_rate = new_water_rate + new_oil_rate + new_gas_rate
    new_control = replace_target(ctrl, ReservoirVoidageTarget(new_rate, new_weights))
    return (new_control, true)
end

