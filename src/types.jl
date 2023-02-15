
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
        new{typeof(equation_of_state), T, O, typeof(reference_densities)}(phases, c, equation_of_state, reference_densities)
    end
end
const LVCompositionalModel = SimulationModel{D, S, F, C} where {D, S<:MultiPhaseCompositionalSystemLV{<:Any, <:Any, <:Any}, F, C}

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
    if has_water
        phase_ind[1] = findfirst(isequal(AqueousPhase()), phases)
        offset = 1
    else
        offset = 0
    end
    phase_ind[1 + offset] = findfirst(isequal(LiquidPhase()), phases)
    phase_ind[2 + offset] = findfirst(isequal(VaporPhase()), phases)
    phase_ind = tuple(phase_ind...)
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

function Base.show(io::IO, d::StandardBlackOilSystem)
    print(io, "StandardBlackOilSystem with $(d.phases)")
end

const BlackOilVariableSwitchingSystem = StandardBlackOilSystem{<:Any, <:Any, <:Any, <:Any, :varswitch, <:Any, <:Any}
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
    function ImmiscibleSystem(phases; reference_densities = ones(length(phases)))
        phases = tuple(phases...)
        reference_densities = tuple(reference_densities...)
        new{typeof(phases), typeof(reference_densities)}(phases, reference_densities)
    end
end
Base.show(io::IO, t::ImmiscibleSystem) = print(io, "ImmiscibleSystem with $(join([typeof(p) for p in t.phases], ", "))")


struct SinglePhaseSystem{P, F} <: MultiPhaseSystem where {P, F<:AbstractFloat}
    phase::P
    rho_ref::F
    function SinglePhaseSystem(phase = LiquidPhase(); reference_density = 1.0)
        return new{typeof(phase), typeof(reference_density)}(phase, reference_density)
    end
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
    crit_ix = findfirst(x -> x > 0, k) - 1
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

abstract type PorousMediumGrid <: AbstractJutulMesh end
abstract type ReservoirGrid <: PorousMediumGrid end

abstract type WellGrid <: PorousMediumGrid
    # Wells are not porous themselves per se, but they are discretizing
    # part of a porous medium.
end

function Base.show(io::IO, w::WellGrid)
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
struct SimpleWell{SC, P} <: WellGrid where {SC, P}
    perforations::P
    surface::SC
    name::Symbol
end

function SimpleWell(reservoir_cells; name = :Well, surface_conditions = default_surface_cond(), kwarg...)
    nr = length(reservoir_cells)
    WI, gdz = common_well_setup(nr; kwarg...)
    perf = (self = ones(Int64, nr), reservoir = vec(reservoir_cells), WI = WI, gdz = gdz)
    return SimpleWell(perf, surface_conditions, name)
end
struct MultiSegmentWell{V, P, N, A, C, SC, S} <: WellGrid
    volumes::V          # One per cell
    perforations::P     # (self -> local cells, reservoir -> reservoir cells, WI -> connection factor)
    neighborship::N     # Well cell connectivity
    top::A              # "Top" node where scalar well quantities live
    centers::C          # Coordinate centers of nodes
    surface::SC         # p, T at surface
    name::Symbol        # Symbol that names the well
    segment_models::S   # Segment pressure drop model for each segment
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
        new{typeof(volumes), typeof(perf), typeof(N), typeof(accumulator), typeof(ext_centers), typeof(surface_conditions), typeof(segment_models)}(volumes, perf, N, accumulator, ext_centers, surface_conditions, name, segment_models)
    end
end

