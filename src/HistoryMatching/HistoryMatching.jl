module HistoryMatching
    using Jutul, JutulDarcy, Jutul.PrettyTables
    export match_injectors!, match_producers!, history_match_objective, evaluate_match
    include("types.jl")
    include("api.jl")
    include("utils.jl")
    include("calculation.jl")
    include("calculation_sum.jl")
end
