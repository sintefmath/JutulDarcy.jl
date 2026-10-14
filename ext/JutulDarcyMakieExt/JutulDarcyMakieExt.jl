module JutulDarcyMakieExt
    using JutulDarcy, Jutul, Makie, Printf, Dates, LinearAlgebra
    import JutulDarcy: well_unit_conversion

    include("well_plots.jl")
    include("multisegment_well_plots.jl")
    include("co2_plots.jl")
    include("faults.jl")
    include("summary.jl")
    include("historymatch.jl")
end
