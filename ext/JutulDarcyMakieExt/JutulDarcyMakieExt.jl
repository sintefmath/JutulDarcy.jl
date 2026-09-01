module JutulDarcyMakieExt
    using JutulDarcy, Jutul, Makie, Printf
    import JutulDarcy: well_unit_conversion
    include("co2_plots.jl")
    include("faults.jl")
    include("summary.jl")
    include("historymatch.jl")
end
