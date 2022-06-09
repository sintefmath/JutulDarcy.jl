
abstract type MultiPhaseSystem <: JutulSystem end
abstract type MultiComponentSystem <: MultiPhaseSystem end
abstract type CompositionalSystem <: MultiComponentSystem end
abstract type BlackOilSystem <: MultiComponentSystem end

abstract type PhaseVariables <: GroupedVariables end
abstract type ComponentVariable <: GroupedVariables end

struct MultiPhaseCompositionalSystemLV{E, T, O} <: CompositionalSystem where T<:Tuple
    phases::T
    components
    equation_of_state::E
    function MultiPhaseCompositionalSystemLV(equation_of_state, phases = (LiquidPhase(), VaporPhase()); other_name = "Water")
        c = copy(equation_of_state.mixture.component_names)
        phases = tuple(phases...)
        T = typeof(phases)
        nph = length(phases)
        @assert nph == 2 || nph == 3
        if nph == 3
            other = only(filter(x -> !(isa(x, LiquidPhase) || isa(x, VaporPhase)), phases))
            O = typeof(other)
            push!(c, other_name)
        else
            O = Nothing
        end
        only(findall(isequal(LiquidPhase()), phases))
        only(findall(isequal(VaporPhase()), phases))
        new{typeof(equation_of_state), T, O}(phases, c, equation_of_state)
    end
end

struct StandardBlackOilSystem{D, W, R, F} <: BlackOilSystem where R<:Real
    saturation_table::D
    rhoLS::R
    rhoVS::R
    function StandardBlackOilSystem(sat; water = true, rhoLS = 1.0, rhoVS = 1.0, formulation::Symbol = :varswitch)
        formulation = :zg
        @assert formulation == :varswitch || formulation == :zg
        new{typeof(sat), water, typeof(rhoLS), formulation}(sat, rhoLS, rhoVS)
    end
end

const BlackOilModelVariableSwitching = SimulationModel{<:Any, StandardBlackOilSystem{D, W, R, :varswitch}, <:Any, <:Any} where {D, W, R}
const BlackOilModelGasFraction = SimulationModel{<:Any, StandardBlackOilSystem{D, W, R, :zg}, <:Any, <:Any} where {D, W, R}
const StandardBlackOilModel = SimulationModel{<:Any, <:StandardBlackOilSystem, <:Any, <:Any}

struct ImmiscibleSystem{T} <: MultiPhaseSystem where T<:Tuple
    phases::T
end
Base.show(io::IO, t::ImmiscibleSystem) = print(io, "ImmiscibleSystem with $(join([typeof(p) for p in t.phases], ", "))")


ImmiscibleSystem(phases::AbstractVector) = ImmiscibleSystem(tuple(phases...))

struct SinglePhaseSystem <: MultiPhaseSystem
    phase
    SinglePhaseSystem(phase = LiquidPhase()) = new(phase)
end
