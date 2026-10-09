"""
Embedded fracture models on an unchanged three-dimensional matrix mesh.
"""
module EFM
    using Jutul, JutulDarcy, LinearAlgebra, StaticArrays
    import JutulDarcy: FractureMatrixConnection

    export AbstractEmbeddedFracture, PolygonFracture, RectangleFracture
    export EllipticalFracture, EllipsoidFracture, fracture_polygon
    export embedded_fracture_domain, embedded_fracture_mesh
    export combined_fracture_domain, setup_embedded_fracture_model
    export EmbeddedFractureMesh, EmbeddedReservoirMesh

    include("fractures.jl")
    include("geometry.jl")
    include("spatial.jl")
    include("mesh.jl")
    include("preprocessing.jl")
    include("domains.jl")
    include("models.jl")
end
