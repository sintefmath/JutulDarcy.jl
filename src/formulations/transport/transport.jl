struct TransportFormulation <: Jutul.JutulFormulation

end

const TransportModel = SimulationModel{<:Any, <:Any, TransportFormulation, <:Any}

include("variables.jl")
include("overloads.jl")
include("functions.jl")
