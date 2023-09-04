
abstract type MultiPhaseSystem <: JutulSystem end
abstract type MultiComponentSystem <: MultiPhaseSystem end
const DarcyFlowModel = SimulationModel{<:Any, <:MultiPhaseSystem, <:Any, <:Any}

abstract type CompositionalSystem <: MultiComponentSystem end
const CompositionalModel = SimulationModel{D, S, F, C} where {D, S<:CompositionalSystem, F, C}

abstract type BlackOilSystem <: MultiComponentSystem end

abstract type PhaseVariables <: VectorVariables end
abstract type ComponentVariable <: VectorVariables end

struct MultiPhaseCompositionalSystemLV{E, T, O, R} <: CompositionalSystem where T<:Tuple
    phases::T
    components::Vector{String}
    equation_of_state::E
    rho_ref::R
end
const LVCompositionalModel = SimulationModel{D, S, F, C} where {D, S<:MultiPhaseCompositionalSystemLV{<:Any, <:Any, <:Any}, F, C}

"""
    MultiPhaseCompositionalSystemLV(equation_of_state, phases = (LiquidPhase(), VaporPhase()); reference_densities = ones(length(phases)), other_name = "Water")

Set up a compositional system for a given `equation_of_state` from `MultiComponentFlash`.
"""
function MultiPhaseCompositionalSystemLV(equation_of_state, phases = (LiquidPhase(), VaporPhase()); reference_densities = ones(length(phases)), other_name = "Water")
    c = copy(equation_of_state.mixture.component_names)
    phases = tuple(phases...)
    T = typeof(phases)
    nph = length(phases)
    @assert nph == 2 || nph == 3
    reference_densities = tuple(reference_densities...)
    @assert length(reference_densities) == nph
    if nph == 3
        other = only(filter(x -> !(isa(x, LiquidPhase) || isa(x, VaporPhase)), phases))
        O = typeof(other)
        push!(c, other_name)
    else
        O = Nothing
    end
    only(findall(isequal(LiquidPhase()), phases))
    only(findall(isequal(VaporPhase()), phases))
    MultiPhaseCompositionalSystemLV{typeof(equation_of_state), T, O, typeof(reference_densities)}(phases, c, equation_of_state, reference_densities)
end


export StandardBlackOilSystem
struct StandardBlackOilSystem{D, V, W, R, F, T, P, Num} <: BlackOilSystem
    rs_max::D
    rv_max::V
    rho_ref::R
    phase_indices::T
    phases::P
    saturated_chop::Bool
    keep_bubble_flag::Bool
    rs_eps::Num
    rv_eps::Num
    s_eps::Num
end

"""
    StandardBlackOilSystem(; rs_max = nothing,
                             rv_max = nothing,
                             phases = (AqueousPhase(), LiquidPhase(), VaporPhase()),
                             reference_densities = [786.507, 1037.84, 0.969758])

Set up a standard black-oil system. Keyword arguments `rs_max` and `rv_max` can
either be nothing or callable objects / functions for the maximum Rs and Rv as a
function of pressure. `phases` can be specified together with
`reference_densities` for each phase. 

NOTE: For the black-oil model, the reference densities significantly impact many
aspects of the PVT behavior. These should generally be set consistently with the
other properties.
"""
function StandardBlackOilSystem(; rs_max::RS = nothing,
                                  rv_max::RV = nothing,
                                  phases = (AqueousPhase(), LiquidPhase(), VaporPhase()),
                                  reference_densities = [786.507, 1037.84, 0.969758], 
                                  saturated_chop = false,
                                  keep_bubble_flag = true,
                                  eps_s = 1e-5,
                                  eps_rs = nothing,
                                  eps_rv = nothing,
                                  formulation::Symbol = :varswitch) where {RS, RV}
    phases = tuple(phases...)
    nph = length(phases)
    if nph == 2 && length(reference_densities) == 3
        reference_densities = reference_densities[2:3]
    end
    reference_densities = tuple(reference_densities...)
    @assert LiquidPhase() in phases
    @assert VaporPhase() in phases
    @assert nph == 2 || nph == 3
    @assert length(reference_densities) == nph
    phase_ind = zeros(Int64, nph)
    has_water = nph == 3
    phase_ind = generate_phase_indices(phases)
    if isnothing(eps_rs)
        if isnothing(rs_max)
            eps_rs = eps_s
        else
            eps_rs = 1e-4*mean(diff(rs_max.F))
        end
    end
    if isnothing(eps_rv)
        if isnothing(rv_max)
            eps_rv = eps_s
        else
            eps_rv = 1e-4*mean(diff(rv_max.F))
        end
    end
    @assert formulation == :varswitch || formulation == :zg
    return StandardBlackOilSystem{RS, RV, has_water, typeof(reference_densities), formulation, typeof(phase_ind), typeof(phases), Float64}(rs_max, rv_max, reference_densities, phase_ind, phases, saturated_chop, keep_bubble_flag, eps_rs, eps_rv, eps_s)
end

function has_other_phase(sys::StandardBlackOilSystem{A, B, W}) where {A, B, W}
    return W
end

function Base.show(io::IO, d::StandardBlackOilSystem)
    print(io, "StandardBlackOilSystem with $(d.phases)")
end

const BlackOilVariableSwitchingSystem = StandardBlackOilSystem{<:Any, <:Any, <:Any, <:Any, :varswitch, <:Any, <:Any}
const BlackOilVariableSwitchingSystemWithWater = StandardBlackOilSystem{<:Any, <:Any, true, <:Any, :varswitch, <:Any, <:Any}
const BlackOilVariableSwitchingSystemWithoutWater = StandardBlackOilSystem{<:Any, <:Any, false, <:Any, :varswitch, <:Any, <:Any}

const BlackOilGasFractionSystem = StandardBlackOilSystem{<:Any, <:Any, <:Any, <:Any, :zg, <:Any, <:Any}

const DisgasBlackOilSystem = StandardBlackOilSystem{<:Any, Nothing, <:Any, <:Any, <:Any, <:Any, <:Any}
const VapoilBlackOilSystem = StandardBlackOilSystem{Nothing, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any}

const BlackOilModelVariableSwitching = SimulationModel{<:Any, BlackOilVariableSwitchingSystem, <:Any, <:Any}
const BlackOilModelGasFraction       = SimulationModel{<:Any, BlackOilGasFractionSystem,        <:Any, <:Any}
const StandardBlackOilModel          = SimulationModel{<:Any, <:StandardBlackOilSystem, <:Any, <:Any}
const VapoilBlackOilModel            = SimulationModel{<:Any, <:VapoilBlackOilSystem, <:Any, <:Any}
const DisgasBlackOilModel            = SimulationModel{<:Any, <:DisgasBlackOilSystem, <:Any, <:Any}

const StandardBlackOilModelWithWater = SimulationModel{<:Any, <:StandardBlackOilSystem{<:Any, <:Any, true, <:Any, <:Any, <:Any, <:Any}, <:Any, <:Any}

struct ImmiscibleSystem{T, F} <: MultiPhaseSystem where {T<:Tuple, F<:NTuple}
    phases::T
    rho_ref::F
end

"""
    ImmiscibleSystem(phases; reference_densities = ones(length(phases)))
    ImmiscibleSystem((LiquidPhase(), VaporPhase()), (1000.0, 700.0))

Set up an immiscible system for the given phases with optional reference
densitites. This system is easy to specify with [Pressure](@ref) and
[Saturations](@ref) as the default primary variables. Immiscible system assume
that there is no mass transfer between phases and that a phase is uniform in
composition.
"""
function ImmiscibleSystem(phases; reference_densities = ones(length(phases)))
    phases = tuple(phases...)
    reference_densities = tuple(reference_densities...)
    return ImmiscibleSystem(phases, reference_densities)
end

Base.show(io::IO, t::ImmiscibleSystem) = print(io, "ImmiscibleSystem with $(join([typeof(p) for p in t.phases], ", "))")


struct SinglePhaseSystem{P, F} <: MultiPhaseSystem where {P, F<:AbstractFloat}
    phase::P
    rho_ref::F
end

"""
    SinglePhaseSystem(phase = LiquidPhase(); reference_density = 1.0)

A single-phase system that only solves for pressure.
"""
function SinglePhaseSystem(phase = LiquidPhase(); reference_density = 1.0)
    return SinglePhaseSystem{typeof(phase), typeof(reference_density)}(phase, reference_density)
end

number_of_components(sys::SinglePhaseSystem) = 1


struct PhaseRelPerm{T, N}
    k::T
    label::Symbol
    connate::N
    critical::N
    s_max::N
    k_max::N
end

export PhaseRelPerm

function PhaseRelPerm(s, k; label = :w, connate = s[1], epsilon = 1e-16)
    for i in eachindex(s)
        if i == 1
            if s[1] == 0.0
                @assert k[i] == 0.0 "Rel. perm. must be zero for s = 0, was $(k[i])"
            end
        else
            msg = "k = $(k[i]) at entry $i corresponding to saturation $(s[i])"
            @assert s[i] >= s[i-1] "Saturations must be increasing: $msg"
            @assert k[i] >= k[i-1] "Rel. Perm. function must be increasing: $msg"
            if k[i] > 1.0
                @warn "Rel. Perm. $label has value larger than 1.0: $msg"
            end
        end
    end
    k_max, ix = findmax(k)
    s_max = s[ix]
    # Last immobile point in table
    crit_ix = findfirst(x -> x > eps(Float64), k) - 1
    crit = s[crit_ix]
    s, k = JutulDarcy.add_missing_endpoints(s, k)
    JutulDarcy.ensure_endpoints!(s, k, epsilon)
    kr = get_1d_interpolator(s, k, cap_endpoints = false)
    return PhaseRelPerm(kr, label, connate, crit, s_max, k_max)
end

(kr::PhaseRelPerm)(S) = kr.k(S)

function Base.show(io::IO, t::MIME"text/plain", kr::PhaseRelPerm)
    println(io, "PhaseRelPerm for $(kr.label):")
    println(io, "  .k: Internal representation: $(kr.k)")
    println(io, "  Connate saturation = $(kr.connate)")
    println(io, "  Critical saturation = $(kr.critical)")
    println(io, "  Maximum rel. perm = $(kr.k_max) at $(kr.s_max)")
end


@enum FlowSourceType begin
    MassSource
    StandardVolumeSource
    VolumeSource
end

struct SourceTerm{I, F, T} <: JutulForce
    cell::I
    value::F
    fractional_flow::T
    type::FlowSourceType
end

export FlowBoundaryCondition
struct FlowBoundaryCondition{I, F, T} <: JutulForce
    cell::I
    pressure::F
    temperature::F
    trans_flow::F
    trans_thermal::F
    fractional_flow::T
    density::Union{F, Nothing}
end

abstract type PorousMediumDomain <: JutulMesh end
abstract type ReservoirGrid <: PorousMediumDomain end

abstract type WellDomain <: PorousMediumDomain
    # Wells are not porous themselves per se, but they are discretizing
    # part of a porous medium.
end

function Base.show(io::IO, w::WellDomain)
    if w isa SimpleWell
        nseg = 0
        nn = 1
        n = "SimpleWell"
    else
        nseg = size(w.neighborship, 2)
        nn = length(w.volumes)
        n = "MultiSegmentWell"
    end
    print(io, "$n [$(w.name)] ($(nn) nodes, $(nseg) segments, $(length(w.perforations.reservoir)) perforations)")
end
struct SimpleWell{SC, P, V} <: WellDomain where {SC, P}
    volume::V
    perforations::P
    surface::SC
    name::Symbol
    explicit_dp::Bool
end

"""
    SimpleWell(reservoir_cells)

Set up a simple well.

NOTE: `setup_vertical_well` or `setup_well` are the recommended way of setting
up wells.
"""
function SimpleWell(
    reservoir_cells;
    name = :Well,
    explicit_dp = true,
    surface_conditions = default_surface_cond(),
    volume = 1000.0, # Regularization volume for well, not a real volume
    kwarg...
    )
    nr = length(reservoir_cells)
    WI, gdz = common_well_setup(nr; kwarg...)
    perf = (self = ones(Int64, nr), reservoir = vec(reservoir_cells), WI = WI, gdz = gdz)
    return SimpleWell(volume, perf, surface_conditions, name, explicit_dp)
end
struct MultiSegmentWell{V, P, N, A, C, SC, S} <: WellDomain
    volumes::V          # One per cell
    perforations::P     # (self -> local cells, reservoir -> reservoir cells, WI -> connection factor)
    neighborship::N     # Well cell connectivity
    top::A              # "Top" node where scalar well quantities live
    centers::C          # Coordinate centers of nodes
    surface::SC         # p, T at surface
    name::Symbol        # Symbol that names the well
    segment_models::S   # Segment pressure drop model for each segment
end

"""
    MultiSegmentWell(reservoir_cells, volumes, centers;
                    N = nothing,
                    name = :Well,
                    perforation_cells = nothing,
                    segment_models = nothing,
                    reference_depth = 0,
                    dz = nothing,
                    surface_conditions = default_surface_cond(),
                    accumulator_volume = mean(volumes),
                    )

Create well perforated in a vector of `reservoir_cells` with corresponding
`volumes` and cell `centers`.

NOTE: `setup_vertical_well` or `setup_well` are the recommended way of setting
up wells.
"""
function MultiSegmentWell(reservoir_cells, volumes::AbstractVector, centers;
                                                    N = nothing,
                                                    name = :Well,
                                                    perforation_cells = nothing,
                                                    segment_models = nothing,
                                                    reference_depth = 0,
                                                    dz = nothing,
                                                    surface_conditions = default_surface_cond(),
                                                    accumulator_volume = mean(volumes),
                                                    kwarg...)
    nv = length(volumes)
    nc = nv + 1
    reservoir_cells = vec(reservoir_cells)
    nr = length(reservoir_cells)
    if isnothing(N)
        @debug "No connectivity. Assuming nicely ordered linear well."
        N = vcat((1:nv)', (2:nc)')
    elseif maximum(N) == nv
        N = vcat([1, 2], N+1)
    end
    nseg = size(N, 2)
    @assert size(N, 1) == 2
    @assert size(centers, 1) == 3
    volumes = vcat([accumulator_volume], volumes)
    ext_centers = hcat([centers[1:2, 1]..., reference_depth], centers)
    @assert length(volumes) == size(ext_centers, 2)
    if !isnothing(reservoir_cells) && isnothing(perforation_cells)
        @assert length(reservoir_cells) == nv "If no perforation cells are given, we must 1->1 correspondence between well volumes and reservoir cells."
        perforation_cells = collect(2:nc)
    end
    perforation_cells = vec(perforation_cells)

    if isnothing(segment_models)
        Δp = SegmentWellBoreFrictionHB(1.0, 1e-4, 0.1)
        segment_models = repeat([Δp], nseg)
    else
        segment_models::AbstractVector
        @assert length(segment_models) == nseg
    end
    if isnothing(dz)
        dz = centers[3, :] - reference_depth
    end
    @assert length(perforation_cells) == nr
    WI, gdz = common_well_setup(nr; dz = dz, kwarg...)
    perf = (self = perforation_cells, reservoir = reservoir_cells, WI = WI, gdz = gdz)
    accumulator = (reference_depth = reference_depth, )
    MultiSegmentWell{typeof(volumes), typeof(perf), typeof(N), typeof(accumulator), typeof(ext_centers), typeof(surface_conditions), typeof(segment_models)}(volumes, perf, N, accumulator, ext_centers, surface_conditions, name, segment_models)
end


export ReservoirSimResult
struct ReservoirSimResult
    wells::AbstractDict
    states::AbstractVector
    time::AbstractVector
    result::Jutul.SimResult
    extra::AbstractDict
end

function ReservoirSimResult(model, result::Jutul.SimResult, forces, extra = Dict(); kwarg...)
    for (k, v) in kwarg
        extra[k] = v
    end
    states, reports = result
    res_states = map(x -> x[:Reservoir], states)
    wells = full_well_outputs(model, states, forces)
    report_time = Jutul.report_times(reports)
    return ReservoirSimResult(wells, res_states, report_time, result, extra)
end

struct TopConditions{N, R}
    density::SVector{N, R}
    volume_fractions::SVector{N, R}
    function TopConditions(n::Int, R::DataType = Float64; density = missing, volume_fractions = missing)
        function internal_convert(x::Missing)
            x = @SVector ones(R, n)
            return x./n
        end
        function internal_convert(x)
            x0 = x
            @assert length(x0) == n
            x = @MVector zeros(R, n)
            @. x = x0
            return SVector(x)
        end
        density = internal_convert(density)
        volume_fractions = internal_convert(volume_fractions)
        return new{n, R}(density, volume_fractions)
    end
end


struct SurfaceWellConditions{T, R} <: ScalarVariable
    storage::T
    separator_conditions::Vector{NamedTuple{(:p, :T), Tuple{R, R}}}
    separator_targets::Vector{Tuple{Int, Int}}
    function SurfaceWellConditions(S::T, c, t, R::DataType = Float64) where T
        new{T, R}(S, c, t)
    end
end

function SurfaceWellConditions(sys::JutulSystem; kwarg...)
    s = Dict{Type, Any}()
    S_t = typeof(default_surface_cond())
    cond = S_t[]
    targets = Tuple{Int, Int}[]
    return SurfaceWellConditions(s, cond, targets)
end
