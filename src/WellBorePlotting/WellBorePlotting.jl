module WellBorePlotting
    using Jutul, JutulDarcy
    import JutulDarcy: MultiSegmentWell, SimpleWell

    include("tubes.jl")
    include("triangulate.jl")
end
