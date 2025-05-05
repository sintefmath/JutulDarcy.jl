struct PressureFormulation <: Jutul.JutulFormulation

end

const PressureModel = SimulationModel{<:Any, <:Any, PressureFormulation, <:Any}

struct PressureReductionFactors <: JutulDarcy.PhaseVariables end

struct PressureEquation{T} <: JutulEquation where T<:ConservationLaw
    conservation::T
end

Jutul.discretization(eq::PressureEquation) = Jutul.discretization(eq.conservation)

include("variables.jl")
include("overloads.jl")
include("functions.jl")
include("multimodel/multimodel_pressure.jl")
