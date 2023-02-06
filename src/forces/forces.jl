
function setup_forces(model::SimulationModel{G, S}; sources = nothing, bc = nothing) where {G<:Any, S<:MultiPhaseSystem}
    if sources isa SourceTerm
        sources = [sources]
    end
    if bc isa FlowBoundaryCondition
        bc = [bc]
    end
    return (sources = sources, bc = bc)
end

include("sources.jl")
include("bc.jl")
