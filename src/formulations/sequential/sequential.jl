struct SequentialSimulator{P, T, S} <: Jutul.JutulSimulator
    pressure::P
    transport::T
    storage::S
end

include("interface.jl")
