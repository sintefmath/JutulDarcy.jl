struct TransportFormulation <: Jutul.JutulFormulation
    schema::Symbol
    function TransportFormulation(; schema::Symbol = :ppu)
        schema in (:ppu, :ppu_nopc, :hybrid) || throw(ArgumentError("Unknown transport schema $schema"))
        new(schema)
    end
end

const TransportModel = SimulationModel{<:Any, <:Any, TransportFormulation, <:Any}

abstract type SequentialFlux <: Jutul.FluxType end

struct TotalSaturationFlux{S} <: SequentialFlux
    function TotalSaturationFlux(schema::Symbol = :ppu)
        schema in (:ppu, :ppu_nopc, :hybrid) || throw(ArgumentError("Unknown flux schema $schema"))
        new{schema}()
    end
end

function TotalSaturationFlux(t::TransportFormulation)
    TotalSaturationFlux(t.schema)
end

transport_scheme(::TotalSaturationFlux{S}) where S = S

include("variables.jl")
include("overloads.jl")
include("upwind.jl")
