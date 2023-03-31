struct PressureFormulation <: Jutul.JutulFormulation

end

const PressureModel = SimulationModel{<:Any, <:Any, PressureFormulation, <:Any}

struct PressureReductionFactors <: JutulDarcy.PhaseVariables end

struct PressureEquation{T} <: JutulEquation where T<:ConservationLaw
    conservation::T
end


include("interface.jl")
include("variables.jl")
include("overloads.jl")
include("functions.jl")