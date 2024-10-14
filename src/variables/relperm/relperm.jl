abstract type AbstractRelativePermeabilities <: PhaseVariables end

function Jutul.default_value(model, v::AbstractRelativePermeabilities)
    # A bit junk, but at least it's junk that sums to one for each cell.
    return 1.0/number_of_phases(model.system)
end

struct ParametricLETRelativePermeabilities <: AbstractRelativePermeabilities

end

struct ParametricCoreyRelativePermeabilities <: AbstractRelativePermeabilities

end

include("hysteresis.jl")
include("endscale.jl")
include("simple.jl")
include("advanced.jl")
