abstract type AbstractRelativePermeabilities <: PhaseVariables end

function Jutul.default_value(model, v::AbstractRelativePermeabilities)
    # A bit junk, but at least it's junk that sums to one for each cell.
    return 1.0/number_of_phases(model.system)
end

Base.@kwdef struct ParametricLETRelativePermeabilities <: AbstractRelativePermeabilities
    wetting_let::Symbol = :WettingLET
    wetting_critical::Symbol = :WettingCritical
    wetting_krmax::Symbol = :WettingKrMax
    nonwetting_let::Symbol = :NonWettingLET
    nonwetting_critical::Symbol = :NonWettingCritical
    nonwetting_krmax::Symbol = :NonWettingKrMax
end

struct ParametricCoreyRelativePermeabilities <: AbstractRelativePermeabilities

end

include("hysteresis.jl")
include("endscale.jl")
include("simple.jl")
include("advanced.jl")
