abstract type AbstractPhase end

"""
Abstract supertype for all multiphase flow systems.
"""
abstract type MultiPhaseSystem <: JutulSystem end
"""
Abstract supertype for multicomponent systems, i.e. flow systems where the
number of components is decoupled from the number of phases.
"""
abstract type MultiComponentSystem <: MultiPhaseSystem end
const DarcyFlowModel = SimulationModel{<:Any, <:MultiPhaseSystem, <:Any, <:Any}

abstract type CompositionalSystem <: MultiComponentSystem end
const CompositionalModel = SimulationModel{D, S, F, C} where {D, S<:CompositionalSystem, F, C}

abstract type BlackOilSystem <: MultiComponentSystem end

abstract type PhaseVariables <: VectorVariables end
abstract type ComponentVariables <: VectorVariables end

struct MultiPhaseCompositionalSystemLV{E, T, O, R, N} <: CompositionalSystem where T<:Tuple
    phases::T
    components::Vector{String}
    equation_of_state::E
    rho_ref::R
    reference_phase_index::Int
end

const LVCompositional2PhaseSystem = MultiPhaseCompositionalSystemLV{<:Any, <:Any, Nothing, <:Any, <:Any}
const LVCompositional3PhaseSystem = MultiPhaseCompositionalSystemLV{<:Any, <:Any, <:AbstractPhase, <:Any, <:Any}

const LVCompositionalModel = SimulationModel{D, S, F, C} where {D, S<:MultiPhaseCompositionalSystemLV{<:Any, <:Any, <:Any, <:Any, <:Any}, F, C}
const LVCompositionalModel2Phase = SimulationModel{D, S, F, C} where {D, S<:LVCompositional2PhaseSystem, F, C}
const LVCompositionalModel3Phase = SimulationModel{D, S, F, C} where {D, S<:LVCompositional3PhaseSystem, F, C}

"""
    MultiPhaseCompositionalSystemLV(equation_of_state)
    MultiPhaseCompositionalSystemLV(equation_of_state, phases = (LiquidPhase(), VaporPhase()); reference_densities = ones(length(phases)), other_name = "Water")

Set up a compositional system for a given `equation_of_state` from
`MultiComponentFlash` with two or three phases. If three phases are provided,
the phase that is not a liquid or a vapor phase will be treated as immiscible in
subsequent simulations and given the name from `other_name` when listed as a
component.
"""
function MultiPhaseCompositionalSystemLV(
        equation_of_state,
        phases = (LiquidPhase(), VaporPhase());
        reference_densities = ones(length(phases)),
        other_name = "Water",
        reference_phase_index = get_reference_phase_index(phases)
    )
    c = copy(MultiComponentFlash.component_names(equation_of_state))
    N = length(c)
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
    return MultiPhaseCompositionalSystemLV{typeof(equation_of_state), T, O, typeof(reference_densities), N}(phases, c, equation_of_state, reference_densities, reference_phase_index)
end

function Base.show(io::IO, sys::MultiPhaseCompositionalSystemLV)
    components = copy(sys.components)
    n = number_of_components(sys)
    if has_other_phase(sys)
        name = "(three-phase)"
        components[end] = "and $(components[end]) as immiscible phase"
        n = n - 1
    else
        name = "(two-phase)"
    end
    eos = sys.equation_of_state
    cnames = join(components, ", ")
    print(io, "MultiPhaseCompositionalSystemLV $name with $(MultiComponentFlash.eostype(eos)) EOS with $n EOS components: $cnames")
end

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
    reference_phase_index::Int
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
function StandardBlackOilSystem(;
        rs_max = nothing,
        rv_max = nothing,
        phases = (AqueousPhase(), LiquidPhase(), VaporPhase()),
        reference_densities = [786.507, 1037.84, 0.969758], 
        saturated_chop = false,
        keep_bubble_flag = true,
        eps_s = 1e-5,
        eps_rs = nothing,
        eps_rv = nothing,
        formulation::Symbol = :varswitch,
        reference_phase_index = missing
    )
    rs_max = region_wrap(rs_max)
    rv_max = region_wrap(rv_max)
    RS = typeof(rs_max)
    RV = typeof(rv_max)
    phases = tuple(phases...)
    if ismissing(reference_phase_index)
        reference_phase_index = get_reference_phase_index(phases)
    end
    reference_phase_index::Int
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
            eps_rs = 1e-4*mean(diff(first(rs_max).F))
        end
    end
    if isnothing(eps_rv)
        if isnothing(rv_max)
            eps_rv = eps_s
        else
            eps_rv = 1e-4*mean(diff(first(rv_max).F))
        end
    end
    @assert formulation == :varswitch || formulation == :zg
    return StandardBlackOilSystem{RS, RV, has_water, typeof(reference_densities), formulation, typeof(phase_ind), typeof(phases), Float64}(rs_max, rv_max, reference_densities, phase_ind, phases, saturated_chop, keep_bubble_flag, eps_rs, eps_rv, eps_s, reference_phase_index)
end

@inline function rs_max_function(sys::StandardBlackOilSystem, region = 1)
    return table_by_region(sys.rs_max, region)
end

@inline function rv_max_function(sys::StandardBlackOilSystem, region = 1)
    return table_by_region(sys.rv_max, region)
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
    reference_phase_index::Int
end

"""
    ImmiscibleSystem(phases; reference_densities = ones(length(phases)))
    ImmiscibleSystem((LiquidPhase(), VaporPhase()), reference_densities = (1000.0, 700.0))

Immiscible flow system: Each component exists only in a single phase, and the
number of components equal the number of phases.

Set up an immiscible system for the given phases with optional reference
densitites. This system is easy to specify with [Pressure](@ref) and
[Saturations](@ref) as the default primary variables. Immiscible system assume
that there is no mass transfer between phases and that a phase is uniform in
composition.
"""
function ImmiscibleSystem(phases; reference_densities = ones(length(phases)), reference_phase_index = missing)
    phases = tuple(phases...)
    if ismissing(reference_phase_index)
        reference_phase_index = get_reference_phase_index(phases)
    end
    reference_densities = tuple(reference_densities...)
    return ImmiscibleSystem(phases, reference_densities, reference_phase_index)
end

Base.show(io::IO, t::ImmiscibleSystem) = print(io, "ImmiscibleSystem with $(join([typeof(p) for p in t.phases], ", "))")


struct SinglePhaseSystem{P, F} <: MultiPhaseSystem where {P, F<:AbstractFloat}
    phase::P
    rho_ref::F
end

const ImmiscibleModel = SimulationModel{D, S, F, C} where {D, S<:ImmiscibleSystem, F, C}

"""
    SinglePhaseSystem(phase = LiquidPhase(); reference_density = 1.0)

A single-phase system that only solves for pressure.
"""
function SinglePhaseSystem(phase = LiquidPhase(); reference_density = 1.0)
    if reference_density isa Real
        reference_density = (reference_density, )
    end
    return SinglePhaseSystem{typeof(phase), typeof(reference_density)}(phase, reference_density)
end

const SinglePhaseModel = SimulationModel{D, S, F, C} where {D, S<:SinglePhaseSystem, F, C}

number_of_components(sys::SinglePhaseSystem) = 1

abstract type AbstractPhaseRelativePermeability{T, N} end

struct PhaseRelativePermeability{T, N} <: AbstractPhaseRelativePermeability{T, N}
    k::T
    label::Symbol
    "Connate saturation"
    connate::N
    "The saturation at which rel. perm. becomes positive"
    critical::N
    "Maximum saturation at which rel. perm. is k_max"
    s_max::N
    "Largest value of rel. perm."
    k_max::N
    "Largest s value in input saturations"
    input_s_max::N
end


"""
    PhaseRelativePermeability(s, k; label = :w, connate = s[1], epsilon = 1e-16)

Type that stores a sorted phase relative permeability table (given as vectors of
equal length `s` and `k`):

``K_r = K(S)``

Optionally, a label for the phase, the connate saturation and a small epsilon
value used to avoid extrapolation can be specified. The return type holds both
the table, the phase context, the autodetected critical and maximum relative
permeability values and can be passed saturation values to evaluate the
underlying function:

```jldoctest
s = range(0, 1, 50)
k = s.^2
kr = PhaseRelativePermeability(s, k)
round(kr(0.5), digits = 2)

# output

0.25
```
"""
function PhaseRelativePermeability(s, k; label = :w, connate = s[1], epsilon = 1e-16)
    s = collect(s)
    k = collect(k)
    msg(i) = "k = $(k[i]) at entry $i corresponding to saturation $(s[i])"
    s, k = saturation_table_handle_defaults(s, k)
    for i in eachindex(s)
        if i == 1
            if s[1] == 0.0
                k[i] == 0.0 || throw(ArgumentError("Rel. perm. must be zero for s = 0, was $(k[i])"))
            end
        else
            s[i] >= s[i-1] || throw(ArgumentError("Saturations must be increasing: $(msg(i))"))
            k[i] >= k[i-1] || throw(ArgumentError("Rel. Perm. function must be increasing: $(msg(i))"))
            if k[i] > 1.0
                @warn "Rel. Perm. $label has value larger than 1.0: $(msg(i))"
            end
        end
    end
    s_max_table = s[end]
    k_max, ix = findmax(k)
    s_max = s[ix]
    # Last immobile point in table
    crit_ix = findfirst(x -> x > eps(Float64), k) - 1
    crit = s[crit_ix]
    s, k = JutulDarcy.add_missing_endpoints(s, k)
    JutulDarcy.ensure_endpoints!(s, k, epsilon)
    kr = get_1d_interpolator(s, k, cap_endpoints = false, constant_dx = false)
    return PhaseRelativePermeability(kr, label, connate, crit, s_max, k_max, s_max_table)
end

(kr::PhaseRelativePermeability)(S) = kr.k(S)

function Base.show(io::IO, t::MIME"text/plain", kr::PhaseRelativePermeability)
    println(io, "PhaseRelativePermeability for $(kr.label):")
    println(io, "  .k: Internal representation: $(kr.k)")
    println(io, "  Connate saturation = $(kr.connate)")
    println(io, "  Critical saturation = $(kr.critical)")
    println(io, "  Maximum rel. perm = $(kr.k_max) at $(kr.s_max)")
end

"""
MassSource: Source is directly interpreted as component masses.
StandardVolumeSource: Source is volume at standard/surface conditions. References densities are used to convert into mass sources.
VolumeSource: Source is volume at in-situ / reservoir conditions.
"""
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

"""
Abstract supertype for all well domains.
"""
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
    reference_depth::V
end

"""
    SimpleWell(reservoir_cells; <keyword arguments>)

Set up a simple well.

# Note

[`setup_vertical_well`](@ref) or [`setup_well`](@ref) are the recommended
way of setting up wells.

# Fields

$FIELDS

"""
function SimpleWell(
        reservoir_cells;
        name = :Well,
        explicit_dp = true,
        surface_conditions = default_surface_cond(),
        reference_depth = 0.0,
        volume = 1000.0, # Regularization volume for well, not a real volume
        kwarg...
    )
    nr = length(reservoir_cells)
    WI, WIth, gdz = common_well_setup(nr; kwarg...)
    T = promote_type(typeof(reference_depth), typeof(volume), eltype(WI), eltype(WIth), eltype(gdz))
    reference_depth = convert(T, reference_depth)
    volume = convert(T, volume)
    WI = T.(WI)
    WIth = T.(WIth)
    gdz = T.(gdz)
    perf = (self = ones(Int64, nr), reservoir = vec(reservoir_cells), WI = WI, WIth = WIth, gdz = gdz)
    return SimpleWell(volume, perf, surface_conditions, name, explicit_dp, reference_depth)
end

struct MultiSegmentWell{V, P, N, A, C, SC, S, M} <: WellDomain
    type::Symbol
    "One of volumes per node (cell)"
    volumes::V
    "(self -> local cells, reservoir -> reservoir cells, WI -> connection factor)"
    perforations::P
    "Well cell connectivity (connections between nodes)"
    neighborship::N
    "Top node where scalar well quantities live"
    top::A
    "End node(s) for the well"
    end_nodes::Vector{Int64}
    "Coordinate centers of nodes"
    centers::C
    "pressure and temperature conditions at surface"
    surface::SC
    "Name of the well as a Symbol"
    name::Symbol
    "Pressure drop model for seg well segment"
    segment_models::S
    "Thermal conductivity of well material"
    material_thermal_conductivity::M
    "Density of well material"
    material_density::M
    "Specific heat capacity of well material"
    material_heat_capacity::M
    "Well void fraction"
    void_fraction::M
end

"""
    MultiSegmentWell(reservoir_cells, volumes, centers;
                    N = nothing,
                    name = :Well,
                    perforation_cells = nothing,
                    segment_models = nothing,
                    segment_length = nothing,
                    reference_depth = 0,
                    dz = nothing,
                    surface_conditions = default_surface_cond(),
                    accumulator_volume = mean(volumes),
                    )

Create well perforated in a vector of `reservoir_cells` with corresponding
`volumes` and cell `centers`.

# Note

[`setup_vertical_well`](@ref) or [`setup_well`](@ref) are the recommended
way of setting up wells.

# Fields

$FIELDS

"""
function MultiSegmentWell(reservoir_cells, volumes::AbstractVector, centers;
            type = :ms,
            accumulator_center = nothing,
            accumulator_volume = mean(volumes),
            N = nothing,
            end_nodes = missing,
            name = :Well,
            perforation_cells = nothing,
            segment_models = nothing,
            reference_depth = nothing,
            dz = nothing,
            segment_length = nothing,
            friction = 1e-4,
            surface_conditions = default_surface_cond(),
            material_thermal_conductivity = 0.0,
            material_heat_capacity = 420.0,
            material_density = 8000.0,
            void_fraction = 1.0,
            extra_perforation_props = NamedTuple(),
            segment_radius = 0.05,
            kwarg...
    )
    if isnothing(reference_depth)
        if isnothing(accumulator_center)
            reference_depth = 0
        else
            reference_depth = accumulator_center[3]
        end
    end
    if isnothing(accumulator_center)
        accumulator_center = [centers[1, 1], centers[2, 1], reference_depth]
    end
    nv = length(volumes)
    nc = nv + 1
    reservoir_cells = vec(reservoir_cells)
    nr = length(reservoir_cells)
    if isnothing(N)
        @debug "No connectivity. Assuming nicely ordered linear well."
        N = vcat((1:nv)', (2:nc)')
    elseif maximum(N) == nv
        N = hcat([1, 2], N.+1)
    end
    if ismissing(end_nodes)
        from_nodes = unique(N[1, :])
        to_nodes = unique(N[2,:])
        end_nodes = setdiff(to_nodes, from_nodes)
    end
    if length(size(centers)) == 3
        @assert size(centers, 3) == 1
        centers = centers[:, :, 1]
    end
    nseg = size(N, 2)
    @assert size(N, 1) == 2
    @assert size(centers, 1) == 3
    volumes = vcat([accumulator_volume], volumes)
    ext_centers = hcat(accumulator_center, centers)
    @assert length(volumes) == size(ext_centers, 2)
    if !isnothing(reservoir_cells) && isnothing(perforation_cells)
        @assert length(reservoir_cells) == nv "If no perforation cells are given, we must 1->1 correspondence between well volumes and reservoir cells."
        perforation_cells = collect(2:nc)
    end
    perforation_cells = vec(perforation_cells)
    if segment_radius isa Real
        segment_radius = fill(segment_radius, nseg)
    end
    @assert length(segment_radius) == nseg "Segment radius must have length equal to number of segments"

    # Process well material properties
    if length(material_thermal_conductivity) == 1
        material_thermal_conductivity = fill(material_thermal_conductivity, nseg)
    end
    @assert length(material_thermal_conductivity) == nseg 
        "Material thermal conductivity must have length equal to number of segments"

    if length(material_heat_capacity) == 1
        material_heat_capacity = fill(material_heat_capacity, nc)
    end
    @assert length(material_heat_capacity) == nc
        "Material heat capacity must have length equal to number of cells (including accumulator node)"

    if length(material_density) == 1
        material_density = fill(material_density, nc)
    end
    @assert length(material_density) == nc
        "Material density must have length equal to number of cells (including accumulator node)"

    if length(void_fraction) == 1
        void_fraction = vcat(1.0, fill(void_fraction, nc-1))
    end
    @assert length(void_fraction) == nc
        "Void fraction must have length equal to number of cells (including accumulator node)"
    @assert void_fraction[1] == 1.0 "Void fraction for accumulator node must be 1.0"

    if isnothing(segment_models)
        segment_models = Vector{SegmentWellBoreFrictionHB{Float64}}()
        for seg in 1:nseg
            l, r = N[:, seg]
            if isnothing(segment_length)
                L = norm(ext_centers[:, l] - ext_centers[:, r], 2)
            else
                if segment_length isa Real
                    L = segment_length
                else
                    L = segment_length[seg]
                end
            end
            diameter = segment_radius[seg]*2
            Δp = SegmentWellBoreFrictionHB(L, friction, diameter)
            push!(segment_models, Δp)
        end
    else
        segment_models::AbstractVector
        @assert length(segment_models) == nseg
    end
    if isnothing(dz)
        dz = centers[3, :] .- reference_depth
    end
    @assert length(perforation_cells) == nr
    WI, WIth, gdz = common_well_setup(nr; dz = dz, kwarg...)
    perf = (self = perforation_cells, reservoir = reservoir_cells, WI = WI, WIth = WIth, gdz = gdz)
    perf = merge(perf, extra_perforation_props)
    for (k, v) in zip(keys(perf), perf)
        @assert length(v) == nr "Perforation property $k must have length equal to number of reservoir cells"
    end
    accumulator = (reference_depth = reference_depth, )
    MultiSegmentWell{typeof(volumes), typeof(perf), typeof(N), typeof(accumulator), typeof(ext_centers), typeof(surface_conditions), typeof(segment_models), typeof(material_thermal_conductivity)}(
        type, volumes, perf, N, accumulator, end_nodes, ext_centers, surface_conditions,
        name, segment_models, material_thermal_conductivity,
        material_density, material_heat_capacity, void_fraction)
end


struct WellResults
    "Vector of time offsets that the reports are given at"
    time::Vector{Float64}
    "Dict-of-dicts that contains the well outputs"
    wells::Dict{Symbol, AbstractDict}
    "Start date for t = 0"
    start_date::Union{Date, Nothing}
    function WellResults(time, well_results, start_date = nothing)
        return new(time, well_results, start_date)
    end
end

Base.pairs(wr::WellResults) = pairs(wr.wells)

struct ReservoirSimResult
    "Well results as a Dict (output from [`full_well_outputs`](@ref))"
    wells::WellResults
    "Reservoir states for each time-step"
    states::AbstractVector
    "The time the states and well solutions are given at"
    time::AbstractVector
    "Raw simulation results with more detailed well results and reports of solution progress"
    result::Jutul.SimResult
    "Dict for holding additional useful data connected to the simulation"
    extra::AbstractDict
end

"""
    ReservoirSimResult(model, result::Jutul.SimResult, forces, extra = Dict(); kwarg...)

Create a specific reservoir simulation results that contains well curves,
reservoir states, and so on. This is the return type from `simulate_reservoir`.

A `ReservoirSimResult` can be unpacked into well solutions, reservoir states and
reporting times:

```julia
res_result::ReservoirSimResult
ws, states, t = res_result
```

# Fields

$FIELDS
"""
function ReservoirSimResult(model, result::Jutul.SimResult, forces, extra = Dict(); start_date = nothing, kwarg...)
    for (k, v) in kwarg
        extra[k] = v
    end
    states, dt, report_ix = Jutul.expand_to_ministeps(result)
    report_time = cumsum(dt)
    res_states = map(x -> x[:Reservoir], states)
    if forces isa Vector
        forces = forces[report_ix]
    end
    wells = full_well_outputs(model, states, forces)
    well_result = WellResults(report_time, wells, start_date)
    return ReservoirSimResult(well_result, res_states, report_time, result, extra)
end

struct TopConditions{N, R}
    density::SVector{N, R}
    volume_fractions::SVector{N, R}
    function TopConditions(density::SVector{N, R}, volume_fractions::SVector{N, R}) where {N, R}
        return new{N, R}(density, volume_fractions)
    end
end

function TopConditions(n::Int, R::DataType = Float64; density = missing, volume_fractions = missing)
    if !ismissing(density)
        R = promote_type(map(typeof, density)..., R)
    end
    if !ismissing(volume_fractions)
        R = promote_type(map(typeof, volume_fractions)..., R)
    end
    return TopConditions(Val(n), Val(R), density, volume_fractions)
end

function TopConditions(::Val{n}, ::Val{R}, density, volume_fractions) where {n, R}
    function internal_convert(x::Missing)
        x = @SVector ones(R, n)
        return x./n
    end
    function internal_convert(x)
        return convert(SVector{n, R}, x)
    end
    density = internal_convert(density)
    volume_fractions = internal_convert(volume_fractions)
    return TopConditions(density, volume_fractions)
end

function JutulDarcy.TopConditions{N, T}(tc::JutulDarcy.TopConditions{N, F}) where {N, T, F}
    density = map(T, tc.density)
    volume_fractions = map(T, tc.volume_fractions)
    return TopConditions(density, volume_fractions)
end

function Base.convert(::Type{TopConditions{N, R}}, tc::TopConditions{N, Float64}) where {N, R}
    return TopConditions(convert.(R, tc.density), convert.(R, tc.volume_fractions))
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

struct PrepareStepWellSolver end


struct PhasePotentials <: PhaseVariables

end

struct AdjustedCellDepths <: ScalarVariable

end

struct CriticalKrPoints <: ScalarVariable end

struct MaxRelPermPoints <: ScalarVariable end

struct LETCoefficients <: JutulVariables end

struct CoreyExponentKrPoints <: ScalarVariable end
