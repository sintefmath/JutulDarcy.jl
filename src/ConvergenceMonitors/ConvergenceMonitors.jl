module ConvergenceMonitors

    using Jutul, JutulDarcy

    export ContractionFactorCuttingCriterion
    export set_contraction_factor_cutting_criterion!

    include("distance_functions.jl")
    include("contraction_factors.jl")
    include("cutting_criterions.jl")
    include("utils.jl")

end