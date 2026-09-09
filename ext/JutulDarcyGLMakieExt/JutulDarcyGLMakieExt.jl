module JutulDarcyGLMakieExt
    using Jutul, JutulDarcy, GLMakie, LinearAlgebra
    import JutulDarcy: well_unit_conversion
    include("well_plots.jl")
    include("multisegment_well_plots.jl")
end
