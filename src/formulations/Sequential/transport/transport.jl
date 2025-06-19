struct TransportFormulation <: Jutul.JutulFormulation

end

const TransportModel = SimulationModel{<:Any, <:Any, TransportFormulation, <:Any}

abstract type SequentialFlux <: Jutul.FluxType end

struct TotalSaturationFlux <: SequentialFlux end

include("variables.jl")
include("overloads.jl")
include("upwind.jl")
