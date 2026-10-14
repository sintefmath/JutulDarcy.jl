module CSP11
    using Jutul, JutulDarcy, StaticArrays, Artifacts, LazyArtifacts

    include("specification.jl")
    include("geometry.jl")
    include("setup.jl")
    include("reporting.jl")

end # module CSP11
