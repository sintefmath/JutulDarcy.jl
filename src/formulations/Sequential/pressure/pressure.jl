function Jutul.setup_equation_storage(model,
        eq::PressureEquation{ConservationLaw{A, B, C, D}}, storage; extra_sparsity = nothing, kwarg...
        ) where {A, B<:TwoPointPotentialFlowHardCoded, C, D}
    return PressureEquationTPFAStorage(model, eq; kwarg...)
end

Jutul.discretization(eq::PressureEquation) = Jutul.discretization(eq.conservation)

include("variables.jl")
include("overloads.jl")
include("functions.jl")
include("bc_and_sources.jl")
include("multimodel/multimodel_pressure.jl")